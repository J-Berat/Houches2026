module DynamoAnalysis

using CairoMakie
import REPL

export BatchConfig, available_figures, run_batch

const PROJECT_DIRECTORY = normpath(joinpath(@__DIR__, ".."))
const MASTER_NOTEBOOK =
    joinpath(PROJECT_DIRECTORY, "dynamo_diagnostics.jl")
const NAVIGATION_CELL_ID =
    Base.UUID("e734297f-506e-45e1-8cb7-b2ae671893eb")
const EARLY_DEFINITION_CELL_IDS = Base.UUID[
    Base.UUID("34d6b5f4-9c2d-42d5-9034-543aeb8ae151"),
    Base.UUID("8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686"),
    # Pure loader/numerical helper definitions. Pluto can make these methods
    # available before their narrative cell position; sequential Julia needs
    # that ordering to be explicit.
    Base.UUID("a8ef96ab-0ddd-4eb2-a216-b7d96c2a9a08"),
]

include(joinpath(@__DIR__, "FigureRegistry.jl"))
include(joinpath(@__DIR__, "BatchCellIndex.jl"))

Base.@kwdef struct BatchConfig
    data_repository::String
    simulations::Vector{String}
    snapshot::Union{Int,Symbol} = :last
    line_of_sight::String = "z"
    figures::Vector{String}
    output_directory::String
    output_format::String = "png"
end

available_figures() = sort(collect(keys(FIGURE_REGISTRY)))
module_isdefined(workspace, name) =
    Base.invokelatest(isdefined, workspace, name)
module_getfield(workspace, name) =
    Base.invokelatest(getfield, workspace, name)
module_names(workspace) =
    Base.invokelatest(names, workspace; all = true)

function source_checksum(path)
    checksum = UInt64(0xcbf29ce484222325)
    for byte in read(path)
        checksum = (checksum ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    end
    string(checksum; base = 16, pad = 16)
end

function validate_source_index()
    actual = source_checksum(MASTER_NOTEBOOK)
    actual == MASTER_SOURCE_CHECKSUM && return
    error(
        "dynamo_diagnostics.jl changed after the batch dependency index was " *
        "generated. Run: julia --startup-file=no --project=. " *
        "tools/generate_batch_index.jl",
    )
end

"""
Read Pluto cells as ordinary Julia source without importing Pluto.

The notebook remains the unique scientific source. Cell UUIDs and their
topological order come from the generated static dependency index.
"""
function read_notebook_cells(path)
    cells = Dict{Base.UUID,String}()
    current_id = nothing
    buffer = IOBuffer()

    function store_current!()
        isnothing(current_id) && return
        cells[current_id] = String(take!(buffer))
    end

    for line in readlines(path; keep = true)
        startswith(line, "# ╔═╡ Cell order:") && begin
            store_current!()
            current_id = nothing
            break
        end
        matched = match(
            r"^# ╔═╡ ([0-9a-fA-F-]{36})\s*$",
            chomp(line),
        )
        if isnothing(matched)
            isnothing(current_id) || write(buffer, line)
        else
            store_current!()
            current_id = Base.UUID(matched.captures[1])
        end
    end
    isnothing(current_id) || store_current!()
    cells
end

function navigation_code(config::BatchConfig)
    simulations_literal = repr(config.simulations)
    snapshot_literal = repr(config.snapshot)
    los_literal = repr(config.line_of_sight)
    """
    begin
        configured_simulations = $simulations_literal

        function configured_run_label(requested_name)
            requested = lowercase(strip(requested_name))
            exact = findfirst(label -> lowercase(label) == requested, run_labels)
            !isnothing(exact) && return run_labels[exact]

            directory_match = findfirst(run_labels) do label
                directory = RUN_DIRS[label]
                lowercase(basename(directory)) == requested ||
                    lowercase(basename(dirname(directory))) == requested
            end
            !isnothing(directory_match) && return run_labels[directory_match]

            error(
                "Simulation introuvable: " * requested_name *
                ". Simulations disponibles: " * join(run_labels, ", ")
            )
        end

        comparison_run_selection =
            unique(configured_run_label.(configured_simulations))
        isempty(comparison_run_selection) &&
            error("La liste des simulations ne doit pas être vide.")

        selected_run = first(comparison_run_selection)
        requested_snapshot = $snapshot_literal
        selected_snapshot = requested_snapshot === :last ?
            length(run_files[selected_run]) :
            clamp(Int(requested_snapshot), 1, length(run_files[selected_run]))
        los_name = $los_literal
        los_name in ("x", "y", "z") ||
            error("La ligne de visée doit être x, y ou z.")
        nothing
    end
    """
end

function bootstrap_workspace()
    workspace_name = gensym(:DynamoBatchWorkspace)
    Core.eval(
        Main,
        Expr(:module, true, workspace_name, Expr(:block)),
    )
    workspace = Base.invokelatest(getfield, Main, workspace_name)
    Base.include_string(
        REPL.softscope,
        workspace,
        """
        using Markdown
        using InteractiveUtils

        macro bind(definition, element)
            quote
                local widget = \$(esc(element))
                global \$(esc(definition)) =
                    Core.applicable(Base.get, widget) ? Base.get(widget) :
                    try
                        local modules = Base.loaded_modules
                        local package = Base.PkgId(
                            Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                            "AbstractPlutoDingetjes",
                        )
                        modules[package].Bonds.initial_value(widget)
                    catch
                        missing
                    end
                widget
            end
        end
        """,
        MASTER_NOTEBOOK,
    )
    workspace
end

function selected_cell_ids(figures)
    selected = Set{Base.UUID}()
    for name in figures
        union!(selected, FIGURE_CELL_IDS[name])
    end
    selected
end

function execute_cells(config::BatchConfig)
    validate_source_index()
    cells = read_notebook_cells(MASTER_NOTEBOOK)
    selected = selected_cell_ids(config.figures)
    missing_cells = setdiff(selected, Set(keys(cells)))
    isempty(missing_cells) ||
        error("Indexed cells missing from the master notebook: $missing_cells")

    workspace = bootstrap_workspace()
    # Pluto predeclares the globals assigned by a cell before applying soft
    # scope. Do the same for requested figure outputs so assignments inside
    # conditional layout branches bind to the workspace module.
    declarations = join(
        [string(FIGURE_REGISTRY[name], " = nothing") for name in config.figures],
        "\n",
    )
    Base.include_string(workspace, declarations, MASTER_NOTEBOOK)
    execution_sequence = unique(vcat(
        [id for id in EARLY_DEFINITION_CELL_IDS if id in selected],
        CELL_EXECUTION_ORDER,
    ))
    executed = 0
    for cell_id in execution_sequence
        cell_id in selected || continue
        code = cell_id == NAVIGATION_CELL_ID ?
            navigation_code(config) : cells[cell_id]
        try
            Base.include_string(
                REPL.softscope,
                workspace,
                code,
                MASTER_NOTEBOOK,
            )
        catch error_value
            println(stderr, "\nErreur dans la cellule ", cell_id, " :")
            showerror(stderr, error_value, catch_backtrace())
            println(stderr)
            rethrow()
        end
        executed += 1
    end
    workspace, executed
end

function validate(config::BatchConfig)
    isdir(config.data_repository) ||
        error("Répertoire de données introuvable: $(config.data_repository)")
    isempty(config.simulations) &&
        error("La liste des simulations ne doit pas être vide.")
    config.line_of_sight in ("x", "y", "z") ||
        error("La ligne de visée doit être x, y ou z.")
    lowercase(config.output_format) in ("png", "pdf") ||
        error("Le format doit être \"png\" ou \"pdf\".")
    unknown = setdiff(config.figures, available_figures())
    isempty(unknown) ||
        error("Figures inconnues: $(join(unknown, ", "))")
end

function save_figures(workspace, config::BatchConfig)
    mkpath(config.output_directory)
    destinations = String[]
    for figure_name in config.figures
        figure_symbol = FIGURE_REGISTRY[figure_name]
        if !module_isdefined(workspace, figure_symbol)
            defined_figures = sort!(
                String.(filter(
                    name -> startswith(String(name), "fig_"),
                    module_names(workspace),
                )),
            )
            error(
                "La figure $figure_name n'a pas été calculée. Variables de " *
                "figure définies: $(join(defined_figures, ", ")).",
            )
        end
        figure = module_getfield(workspace, figure_symbol)
        isnothing(figure) &&
            error("La figure $figure_name est restée indéfinie après le calcul.")
        destination = joinpath(
            config.output_directory,
            "$(figure_name).$(lowercase(config.output_format))",
        )
        CairoMakie.save(destination, figure)
        push!(destinations, destination)
        println("Figure enregistrée : ", destination)
    end
    destinations
end

function run_batch(config::BatchConfig)
    validate(config)
    if isempty(config.figures)
        println("Figures disponibles :")
        foreach(name -> println("  \"", name, "\","), available_figures())
        return String[]
    end

    ENV["DYNAMO_DATA_REPOSITORY"] = abspath(config.data_repository)
    get!(
        ENV,
        "DYNAMO_RAW_CUBE_CACHE_ENTRIES",
        string(max(1, length(config.simulations))),
    )

    println("Moteur batch : Julia natif (Pluto non chargé)")
    println("Source scientifique unique : ", MASTER_NOTEBOOK)
    println("Simulations : ", join(config.simulations, ", "))
    println("Figures : ", join(config.figures, ", "))
    println(
        "Cache I/O : ",
        ENV["DYNAMO_RAW_CUBE_CACHE_ENTRIES"],
        " cube(s), avec limite mémoire automatique",
    )
    if Threads.nthreads() == 1
        println(
            "Conseil : lancer Julia avec --threads=auto pour paralléliser ",
            "les balayages temporels.",
        )
    else
        println("Threads Julia : ", Threads.nthreads())
    end

    workspace, executed = execute_cells(config)
    println("Cellules scientifiques exécutées : ", executed)
    save_figures(workspace, config)
end

end
