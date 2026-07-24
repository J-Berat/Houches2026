module DynamoAnalysis

using CairoMakie
using Printf
import REPL

export BatchConfig, available_figures, available_notebooks,
    figures_for_notebooks, run_batch

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

function format_duration(seconds)
    seconds < 60 && return @sprintf("%.1f s", seconds)
    minutes, remaining_seconds = divrem(round(Int, seconds), 60)
    minutes < 60 &&
        return @sprintf("%d min %02d s", minutes, remaining_seconds)
    hours, remaining_minutes = divrem(minutes, 60)
    @sprintf("%d h %02d min %02d s", hours, remaining_minutes, remaining_seconds)
end

terminal_is_interactive() =
    lowercase(get(ENV, "TERM", "dumb")) != "dumb" &&
    lowercase(get(ENV, "CI", "false")) ∉ ("1", "true", "yes")

"""
Display a compact progress bar in an interactive terminal.

When stdout is redirected to a log file, emit bounded progress messages instead
of carriage-return animations so the resulting log remains readable.
"""
function show_progress(completed, total, label; force = false)
    total > 0 || return
    completed = clamp(completed, 0, total)
    fraction = completed / total
    percentage = round(Int, 100 * fraction)
    bar_width = 30
    filled = clamp(floor(Int, bar_width * fraction), 0, bar_width)
    bar = repeat("█", filled) * repeat("░", bar_width - filled)

    if terminal_is_interactive()
        print(
            stdout,
            "\r\e[2K[", bar, "] ",
            lpad(percentage, 3), "%  ",
            label,
        )
        completed == total && println(stdout)
        flush(stdout)
        return
    end

    reporting_step = max(1, cld(total, 10))
    if force || completed == 0 || completed == total ||
            completed % reporting_step == 0
        println(
            stdout,
            "[", lpad(percentage, 3), "%] ",
            label,
            " (", completed, "/", total, ")",
        )
    end
end

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
                "Simulation not found: " * requested_name *
                ". Available simulations: " * join(run_labels, ", ")
            )
        end

        comparison_run_selection =
            unique(configured_run_label.(configured_simulations))
        isempty(comparison_run_selection) &&
            error("The simulation list must not be empty.")

        selected_run = first(comparison_run_selection)
        requested_snapshot = $snapshot_literal
        selected_snapshot = requested_snapshot === :last ?
            length(run_files[selected_run]) :
            clamp(Int(requested_snapshot), 1, length(run_files[selected_run]))
        los_name = $los_literal
        los_name in ("x", "y", "z") ||
            error("The line of sight must be x, y, or z.")
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

function cell_description(cell_id, code, config)
    cell_id == NAVIGATION_CELL_ID &&
        return "Selecting simulations and snapshot"
    cell_id in EARLY_DEFINITION_CELL_IDS &&
        return "Initializing scientific functions"

    for figure_name in config.figures
        figure_variable = string(FIGURE_REGISTRY[figure_name])
        assignment = Regex("\\b" * figure_variable * "\\s*=")
        occursin(assignment, code) &&
            return "Building figure \"$figure_name\""
    end

    phases = [
        "comparison_cubes" => "Reading comparison cubes",
        "temporal_run_series" => "Computing time series",
        "power_spectrum" => "Computing power spectra",
        "structure_function" => "Computing structure functions",
        "dust_" => "Computing dust observables",
        "starlight_" => "Computing starlight polarization",
        "zeeman_" => "Computing Zeeman observables",
        "moose_" => "Computing MOOSE observables",
        "shine_" => "Computing SHINE observables",
        "load_raw_cube" => "Reading and converting cubes",
    ]
    for (needle, description) in phases
        occursin(needle, code) && return description
    end
    "Computing scientific dependencies"
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
    planned_cells = [id for id in execution_sequence if id in selected]
    total_cells = length(planned_cells)
    println("Scientific plan: ", total_cells, " required cells")

    executed = 0
    for cell_id in planned_cells
        code = cell_id == NAVIGATION_CELL_ID ?
            navigation_code(config) : cells[cell_id]
        description = cell_description(cell_id, code, config)
        show_progress(executed, total_cells, description; force = true)
        try
            Base.include_string(
                REPL.softscope,
                workspace,
                code,
                MASTER_NOTEBOOK,
            )
        catch error_value
            terminal_is_interactive() && println(stdout)
            println(stderr, "\nError in cell ", cell_id, ":")
            showerror(stderr, error_value, catch_backtrace())
            println(stderr)
            rethrow()
        end
        executed += 1
        show_progress(executed, total_cells, description)
    end
    workspace, executed
end

function validate(config::BatchConfig)
    isdir(config.data_repository) ||
        error("Data directory not found: $(config.data_repository)")
    isempty(config.simulations) &&
        error("The simulation list must not be empty.")
    config.line_of_sight in ("x", "y", "z") ||
        error("The line of sight must be x, y, or z.")
    lowercase(config.output_format) in ("png", "pdf") ||
        error("The output format must be \"png\" or \"pdf\".")
    unknown = setdiff(config.figures, available_figures())
    isempty(unknown) ||
        error("Unknown figures: $(join(unknown, ", "))")
end

function save_figures(workspace, config::BatchConfig)
    output_directory = abspath(config.output_directory)
    mkpath(output_directory)
    destinations = String[]
    total_figures = length(config.figures)
    println("\nSaving ", total_figures, " figure(s)")
    for (index, figure_name) in enumerate(config.figures)
        show_progress(
            index - 1,
            total_figures,
            "Writing \"$figure_name\"";
            force = true,
        )
        figure_symbol = FIGURE_REGISTRY[figure_name]
        if !module_isdefined(workspace, figure_symbol)
            defined_figures = sort!(
                String.(filter(
                    name -> startswith(String(name), "fig_"),
                    module_names(workspace),
                )),
            )
            error(
                "Figure $figure_name was not computed. Defined figure " *
                "variables: $(join(defined_figures, ", ")).",
            )
        end
        figure = module_getfield(workspace, figure_symbol)
        isnothing(figure) &&
            error("Figure $figure_name remained undefined after computation.")
        destination = joinpath(
            output_directory,
            "$(figure_name).$(lowercase(config.output_format))",
        )
        CairoMakie.save(destination, figure)
        push!(destinations, destination)
        show_progress(
            index,
            total_figures,
            "Saved figure: $figure_name",
        )
    end
    destinations
end

function run_batch(config::BatchConfig)
    started_at = time_ns()
    validate(config)
    if isempty(config.figures)
        println("Available figures:")
        foreach(name -> println("  \"", name, "\","), available_figures())
        return String[]
    end

    ENV["DYNAMO_DATA_REPOSITORY"] = abspath(config.data_repository)
    get!(
        ENV,
        "DYNAMO_RAW_CUBE_CACHE_ENTRIES",
        string(max(1, length(config.simulations))),
    )

    output_directory = abspath(config.output_directory)
    println("\n", repeat("═", 72))
    println("DYNAMO — BATCH COMPUTATION")
    println(repeat("═", 72))
    println("Engine              : Native Julia (Pluto not loaded)")
    println("Scientific source   : ", MASTER_NOTEBOOK)
    println("Data directory      : ", ENV["DYNAMO_DATA_REPOSITORY"])
    println("Snapshot            : ", config.snapshot)
    println("Line of sight       : ", config.line_of_sight)
    println("Output format       : ", uppercase(config.output_format))
    println("Output directory    : ", output_directory)
    println("\nCompared simulations (", length(config.simulations), "):")
    foreach(name -> println("  • ", name), config.simulations)
    println("\nRequested figures (", length(config.figures), "):")
    foreach(name -> println("  • ", name), config.figures)
    println(
        "\nI/O cache           : ",
        ENV["DYNAMO_RAW_CUBE_CACHE_ENTRIES"],
        " cube(s), with an automatic memory limit",
    )
    if Threads.nthreads() == 1
        println(
            "Tip: start Julia with --threads=auto to parallelize ",
            "time-series sweeps.",
        )
    else
        println("Julia threads       : ", Threads.nthreads())
    end

    println("\nStarting scientific computations")
    workspace, executed = execute_cells(config)
    destinations = save_figures(workspace, config)
    elapsed_seconds = (time_ns() - started_at) / 1.0e9

    println("\n", repeat("═", 72))
    println("COMPUTATION COMPLETE")
    println(repeat("═", 72))
    println("Executed cells      : ", executed)
    println("Generated figures   : ", length(destinations))
    println("Total runtime       : ", format_duration(elapsed_seconds))
    println("Output directory    : ", output_directory)
    println("\nCreated files:")
    foreach(path -> println("  ✓ ", path), destinations)
    println(repeat("═", 72))
    destinations
end

end
