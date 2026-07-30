import Pkg

const PROJECT_DIRECTORY = @__DIR__
Pkg.activate(PROJECT_DIRECTORY)

include(joinpath(PROJECT_DIRECTORY, "src", "DynamoAnalysis.jl"))
using .DynamoAnalysis

# =============================================================================
# CONFIGURATION TO EDIT
# =============================================================================

const DEFAULT_COMPARISON_REPOSITORY =
    "/Xnfs/Houches2026/DynSim/cooling_freq_output"

# Select any subset of: "mach", "resolution", and "ratio".
const SELECTED_COMPARISONS = [
    "mach",
    "resolution",
    "ratio",
]

# Every projected or synthetic-observation figure is generated for these
# viewing axes. Figures that do not depend on the LOS are generated once in
# the sibling `common` directory.
const SELECTED_LINES_OF_SIGHT = [
    "x",
    "y",
    "z",
]

# A single snapshot count is shared by every selected figure family. Snapshots
# are distributed uniformly across the available run so temporal diagnostics
# cover the complete evolution instead of applying notebook-specific
# first/last rules.
const SELECTED_SNAPSHOT_COUNT = 10

const COMPARISONS = [
    (
        key = "mach",
        folder = "VaryingMach",
        output_name = "varying_mach",
        simulations = [
            "turb_rms_10_N128",
            "turb_rms_50_N128",
            "turb_rms_100_N128",
        ],
    ),
    (
        key = "resolution",
        folder = "VaryingRes",
        output_name = "varying_resolution",
        simulations = [
            "turb_rms_50_N128",
            "turb_rms_50_N256",
        ],
    ),
    (
        key = "ratio",
        folder = "VaryingRatio",
        output_name = "varying_ratio",
        simulations = [
            "turb_rms_50_N128_R0",
            "turb_rms_50_N128_R1",
        ],
    ),
]

# Compute every registered notebook group by default. Replace this expression
# with a string vector to run only a subset. MOOSE deliberately contributes
# only its phi-space tomography and H I--Faraday HOG products.
const SELECTED_NOTEBOOKS = available_notebooks()

const COMPARISON_REPOSITORY = get(
    ENV,
    "DYNAMO_COMPARISON_REPOSITORY",
    DEFAULT_COMPARISON_REPOSITORY,
)

const UNKNOWN_COMPARISONS = setdiff(
    SELECTED_COMPARISONS,
    getproperty.(COMPARISONS, :key),
)
isempty(UNKNOWN_COMPARISONS) || error(
    "Unknown comparisons: $(join(UNKNOWN_COMPARISONS, ", ")). " *
    "Available comparisons: " * join(getproperty.(COMPARISONS, :key), ", "),
)

const UNKNOWN_LINES_OF_SIGHT =
    setdiff(lowercase.(SELECTED_LINES_OF_SIGHT), ["x", "y", "z"])
isempty(UNKNOWN_LINES_OF_SIGHT) || error(
    "Unknown lines of sight: $(join(UNKNOWN_LINES_OF_SIGHT, ", ")). " *
    "Available values: x, y, z.",
)

function batch_config(
        comparison,
        figures,
        snapshot_window,
        snapshot_count,
        output_group,
        output_subdirectory,
        line_of_sight,
        ;
        comparison_repository = COMPARISON_REPOSITORY,
        output_root = joinpath(PROJECT_DIRECTORY, "figures"),
        output_format = "png",
    )
    BatchConfig(
        data_repository = isempty(comparison.folder) ?
            comparison_repository :
            joinpath(comparison_repository, comparison.folder),
        simulations = comparison.simulations,
        snapshot = :last,
        snapshot_window = snapshot_window,
        snapshot_count = snapshot_count,
        line_of_sight = line_of_sight,
        figures = figures,
        output_directory = joinpath(
            output_root,
            comparison.output_name,
            "$(output_group)_$(snapshot_window)$(snapshot_count)",
            output_subdirectory,
        ),
        output_format = output_format,
    )
end

function append_figure_jobs!(
        configurations,
        comparison,
        figures,
        snapshot_window,
        snapshot_count,
        output_group,
        lines_of_sight = SELECTED_LINES_OF_SIGHT;
        comparison_repository = COMPARISON_REPOSITORY,
        output_root = joinpath(PROJECT_DIRECTORY, "figures"),
        output_format = "png",
    )
    split = split_figures_by_los(figures)

    # LOS-independent figures are evaluated exactly once. The nominal z axis
    # only initializes notebook navigation and cannot affect these products.
    isempty(split.independent) || push!(
        configurations,
        batch_config(
            comparison,
            split.independent,
            snapshot_window,
            snapshot_count,
            output_group,
            "common",
            "z",
            comparison_repository = comparison_repository,
            output_root = output_root,
            output_format = output_format,
        ),
    )

    for line_of_sight in lines_of_sight
        isempty(split.dependent) && continue
        push!(
            configurations,
            batch_config(
                comparison,
                split.dependent,
                snapshot_window,
                snapshot_count,
                output_group,
                "los_$(line_of_sight)",
                line_of_sight,
                comparison_repository = comparison_repository,
                output_root = output_root,
                output_format = output_format,
            ),
        )
    end
    configurations
end

const BATCH_CONFIGS = let
    configurations = BatchConfig[]
    for comparison in COMPARISONS
        comparison.key in SELECTED_COMPARISONS || continue
        append_figure_jobs!(
            configurations,
            comparison,
            figures_for_notebooks(SELECTED_NOTEBOOKS),
            :even,
            SELECTED_SNAPSHOT_COUNT,
            "selected",
        )
    end
    configurations
end

# =============================================================================
# EXECUTION — DO NOT EDIT BELOW THIS LINE
# =============================================================================

function prompt_text(question, default)
    print(question, " [", default, "]: ")
    flush(stdout)
    answer = strip(readline(stdin))
    isempty(answer) ? String(default) : answer
end

function selected_indices(answer, count)
    normalized = lowercase(strip(answer))
    normalized in ("", "all", "*") && return collect(1:count)
    selected = Int[]
    for token in split(normalized, ',')
        part = strip(token)
        if occursin('-', part)
            bounds = split(part, '-'; limit = 2)
            length(bounds) == 2 || error("Invalid range: $(part)")
            first_index, last_index = parse.(Int, bounds)
            first_index <= last_index || error("Invalid range: $(part)")
            append!(selected, first_index:last_index)
        else
            push!(selected, parse(Int, part))
        end
    end
    all(index -> 1 <= index <= count, selected) ||
        error("Selection must contain indices between 1 and $(count).")
    unique(selected)
end

function prompt_multiple(title, options; default = "all", labels = options)
    println("\n", title)
    for (index, label) in enumerate(labels)
        println("  ", lpad(index, 2), ". ", label)
    end
    while true
        print("Selection (comma/range/all) [", default, "]: ")
        flush(stdout)
        answer = strip(readline(stdin))
        isempty(answer) && (answer = default)
        try
            indices = selected_indices(answer, length(options))
            isempty(indices) && error("Select at least one entry.")
            return options[indices]
        catch error_value
            println(stderr, "Invalid selection: ", sprint(showerror, error_value))
        end
    end
end

function prompt_one(title, options; default = 1, labels = options)
    println("\n", title)
    for (index, label) in enumerate(labels)
        println("  ", index == default ? " • " : "   ", index, ". ", label)
    end
    while true
        print("Selection [", default, "]: ")
        flush(stdout)
        answer = strip(readline(stdin))
        isempty(answer) && return options[default]
        try
            index = parse(Int, answer)
            1 <= index <= length(options) ||
                error("Choose a number between 1 and $(length(options)).")
            return options[index]
        catch error_value
            println(stderr, "Invalid selection: ", sprint(showerror, error_value))
        end
    end
end

function prompt_positive_integer(question, default)
    while true
        answer = prompt_text(question, string(default))
        try
            value = parse(Int, answer)
            value > 0 || error("The value must be positive.")
            return value
        catch error_value
            println(stderr, "Invalid value: ", sprint(showerror, error_value))
        end
    end
end

function normalized_path(path)
    stripped = strip(path)
    expanded = stripped == "~" ? homedir() :
        startswith(stripped, "~/") ? joinpath(homedir(), stripped[3:end]) :
        stripped
    abspath(expanded)
end

const SNAPSHOT_EXTENSIONS = Set([".h5", ".hdf5", ".fits", ".fit", ".fts"])

function directory_has_direct_snapshots(directory)
    isdir(directory) || return false
    isdir(joinpath(directory, "DataCubes")) && return true
    !isempty(ramses_output_directories(directory)) && return true
    try
        any(readdir(directory; join = true)) do entry
            isfile(entry) && lowercase(splitext(entry)[2]) in SNAPSHOT_EXTENSIONS
        end
    catch
        false
    end
end

function ramses_output_directories(simulation_directory)
    isdir(simulation_directory) || return String[]
    entries = try
        readdir(simulation_directory; join = true, sort = true)
    catch
        return String[]
    end
    sort(filter(entries) do entry
        isdir(entry) || return false
        matched = match(r"^output_(\d+)$", basename(entry))
        isnothing(matched) && return false
        isfile(joinpath(entry, "info_$(matched.captures[1]).txt"))
    end)
end

function immediate_simulation_directories(directory)
    entries = try
        readdir(directory; join = true, sort = true)
    catch error_value
        error(
            "Cannot read $(directory): " *
            sprint(showerror, error_value),
        )
    end
    ignored = Set([
        ".dynamo_cache",
        "figures",
        "exports",
        "powerspectra",
        "power_spectra",
    ])
    filter(entries) do entry
        isdir(entry) &&
            !startswith(basename(entry), ".") &&
            lowercase(basename(entry)) ∉ ignored
    end
end

function safe_output_name(path)
    name = lowercase(basename(normpath(path)))
    cleaned = replace(name, r"[^a-z0-9]+" => "_")
    isempty(strip(cleaned, '_')) ? "custom_comparison" : strip(cleaned, '_')
end

function fits_field_directory_is_snapshot(directory)
    isdir(directory) || return false
    stems = try
        Set(lowercase(splitext(entry)[1]) for entry in readdir(directory)
            if lowercase(splitext(entry)[2]) in (".fits", ".fit", ".fts"))
    catch
        return false
    end
    density_names = Set(["rho", "density", "massdensity"])
    velocity_names = (
        Set(["vx", "velx", "velocityx", "xvelocity"]),
        Set(["vy", "vely", "velocityy", "yvelocity"]),
        Set(["vz", "velz", "velocityz", "zvelocity"]),
    )
    !isempty(intersect(stems, density_names)) &&
        all(!isempty(intersect(stems, names)) for names in velocity_names)
end

function direct_snapshot_sources(directory)
    isdir(directory) || return String[]
    entries = try
        readdir(directory; join = true, sort = true)
    catch error_value
        error(
            "Cannot list snapshots in $(directory): " *
            sprint(showerror, error_value),
        )
    end
    sort(filter(entries) do entry
        (isfile(entry) &&
            lowercase(splitext(entry)[2]) in SNAPSHOT_EXTENSIONS) ||
            fits_field_directory_is_snapshot(entry)
    end)
end

function simulation_snapshot_inventory(
        comparison_repository,
        comparison,
        simulation,
    )
    group_directory = isempty(comparison.folder) ?
        comparison_repository :
        joinpath(comparison_repository, comparison.folder)
    simulation_directory = joinpath(group_directory, simulation)
    isdir(simulation_directory) || error(
        "Simulation directory not found: $(simulation_directory)",
    )
    cube_directory = isdir(joinpath(simulation_directory, "DataCubes")) ?
        joinpath(simulation_directory, "DataCubes") :
        simulation_directory
    sources = direct_snapshot_sources(cube_directory)
    if isempty(sources)
        sources = ramses_output_directories(simulation_directory)
        isempty(sources) || return (
            simulation,
            cube_directory = simulation_directory,
            sources,
            format = "RAMSES",
        )
    end
    isempty(sources) && error(
        "No HDF5 or FITS snapshots were found for simulation $(simulation) " *
        "in $(cube_directory).",
    )
    (; simulation, cube_directory, sources, format = "HDF5/FITS")
end

function print_snapshot_inventory(comparisons, comparison_repository)
    println("\n", repeat("─", 72))
    println("SNAPSHOT INVENTORY")
    println(repeat("─", 72))
    inventories = NamedTuple[]
    for comparison in comparisons
        for simulation in comparison.simulations
            inventory = simulation_snapshot_inventory(
                comparison_repository,
                comparison,
                simulation,
            )
            push!(inventories, inventory)
            println(
                "\n",
                inventory.simulation,
                ": ",
                length(inventory.sources),
                " snapshot(s) [",
                inventory.format,
                "]",
            )
            println("  Cube directory: ", inventory.cube_directory)
            for (index, source) in enumerate(inventory.sources)
                println("  ", lpad(index, 4), ". ", source)
            end
        end
    end
    println(repeat("─", 72))
    inventories
end

function inferred_ramses_resolution(simulation_paths)
    for path in reverse(simulation_paths)
        explicit = match(
            r"(?:^|_)[Nn](\d{2,4})(?=$|_)",
            basename(normpath(path)),
        )
        isnothing(explicit) ||
            return parse(Int, explicit.captures[1])
    end
    resolutions = Int[]
    for path in simulation_paths
        for matched in eachmatch(r"(?:^|_)(\d{2,4})(?=$|_)", basename(normpath(path)))
            push!(resolutions, parse(Int, matched.captures[1]))
        end
    end
    isempty(resolutions) ? 256 : last(resolutions)
end

function ramses_cache_comparison(simulation_paths, requested_path)
    inventories = Dict(
        path => ramses_output_directories(path) for path in simulation_paths)
    all(outputs -> !isempty(outputs), values(inventories)) || error(
        "A custom comparison cannot currently mix raw RAMSES outputs with " *
        "HDF5/FITS simulations. Convert the RAMSES simulations first or " *
        "select simulations using the same format.",
    )

    println("\n", repeat("─", 72))
    println("RAW RAMSES OUTPUTS")
    println(repeat("─", 72))
    selected_outputs = Dict{String,Vector{String}}()
    for path in simulation_paths
        outputs = inventories[path]
        println("\n", basename(normpath(path)), ": ", length(outputs), " output(s)")
        selected_outputs[path] = prompt_multiple(
            "RAMSES outputs to cache for $(basename(normpath(path)))",
            outputs;
            labels = basename.(outputs),
        )
    end

    dry_run = lowercase(get(ENV, "DYNAMO_DRY_RUN", "false")) in
        ("1", "true", "yes")
    if dry_run
        println("\nDry run: RAMSES conversion was not started.")
        comparison = (
            key = "custom",
            folder = "",
            output_name = safe_output_name(requested_path),
            simulations = basename.(normpath.(simulation_paths)),
        )
        return [comparison], dirname(first(simulation_paths))
    end

    resolution = prompt_positive_integer(
        "Uniform RAMSES cache resolution per axis",
        inferred_ramses_resolution(simulation_paths),
    )
    cache_root = normalized_path(prompt_text(
        "RAMSES HDF5 cache root",
        joinpath(PROJECT_DIRECTORY, ".dynamo_cache", "ramses_cubes"),
    ))
    mkpath(cache_root)
    converter = joinpath(PROJECT_DIRECTORY, "tools", "ramses_to_hdf5.py")
    isfile(converter) || error("RAMSES converter not found: $(converter)")
    python = get(ENV, "DYNAMO_PYTHON", "python3")

    cache_names = String[]
    used_names = Set{String}()
    for (simulation_index, path) in enumerate(simulation_paths)
        base_name = safe_output_name(path)
        cache_name = base_name
        suffix = 2
        while cache_name in used_names
            cache_name = "$(base_name)_$(suffix)"
            suffix += 1
        end
        push!(used_names, cache_name)
        push!(cache_names, cache_name)
        cache_directory = joinpath(cache_root, cache_name, "DataCubes")
        mkpath(cache_directory)
        outputs_argument = join(basename.(selected_outputs[path]), ",")
        command = `$(python) $(converter) --simulation $(path) --output-directory $(cache_directory) --resolution $(resolution) --outputs $(outputs_argument)`
        println(
            "\nConverting RAMSES simulation ",
            simulation_index,
            "/",
            length(simulation_paths),
            ": ",
            path,
        )
        try
            run(command)
        catch error_value
            error(
                "RAMSES conversion failed. Install Python dependencies with " *
                "`$(python) -m pip install --user -r " *
                "requirements-ramses.txt`. Original error: " *
                sprint(showerror, error_value),
            )
        end
    end

    comparison = (
        key = "custom",
        folder = "",
        output_name = safe_output_name(requested_path),
        simulations = cache_names,
    )
    [comparison], cache_root
end

function custom_comparison_selection()
    requested_path = normalized_path(prompt_text(
        "Simulation or simulation-group directory",
        COMPARISON_REPOSITORY,
    ))
    if !isdir(requested_path)
        if startswith(requested_path, "/Xnfs/") && Sys.isapple()
            remote_launcher = joinpath(
                PROJECT_DIRECTORY,
                "run_figures_interactive_cluster.sh",
            )
            isfile(remote_launcher) || error(
                "Cluster launcher not found: $(remote_launcher)",
            )
            println(
                "\nThe selected /Xnfs path is not mounted on this laptop.",
            )
            println(
                "Switching automatically to the interactive cluster launcher...",
            )
            flush(stdout)
            withenv("DYNAMO_REMOTE_DATA_PATH" => requested_path) do
                run(`bash $(remote_launcher)`)
            end
            exit(0)
        end
        cluster_hint = startswith(requested_path, "/Xnfs/") ?
            "\nThis is a cluster-only /Xnfs path. Run " *
            "`bash run_figures_interactive_cluster.sh` from the laptop so " *
            "the menu executes through SSH on PSMN_sr650node230." : ""
        error("Data directory not found: $(requested_path)$(cluster_hint)")
    end

    if directory_has_direct_snapshots(requested_path)
        if !isempty(ramses_output_directories(requested_path))
            return ramses_cache_comparison([requested_path], requested_path)
        end
        simulation_path = lowercase(basename(normpath(requested_path))) ==
            "datacubes" ? dirname(normpath(requested_path)) : requested_path
        simulation_name = basename(normpath(simulation_path))
        comparison = (
            key = "custom",
            folder = "",
            output_name = safe_output_name(simulation_path),
            simulations = [simulation_name],
        )
        println("\nDetected one simulation: ", simulation_name)
        return [comparison], dirname(normpath(simulation_path))
    end

    candidates = immediate_simulation_directories(requested_path)
    isempty(candidates) && error(
        "No simulation directory was found directly inside $(requested_path). " *
        "Select the simulation directory itself or its immediate parent.",
    )
    candidate_names = basename.(candidates)
    selected_names = prompt_multiple(
        "Simulations to compare",
        candidate_names;
        labels = candidate_names,
    )
    selected_paths = [joinpath(requested_path, name) for name in selected_names]
    raw_flags = .!isempty.(ramses_output_directories.(selected_paths))
    if any(raw_flags)
        all(raw_flags) || error(
            "The selected simulations mix raw RAMSES and HDF5/FITS data. " *
            "Select only raw RAMSES simulations for automatic conversion.",
        )
        return ramses_cache_comparison(selected_paths, requested_path)
    end
    comparison = (
        key = "custom",
        folder = "",
        output_name = safe_output_name(requested_path),
        simulations = selected_names,
    )
    [comparison], requested_path
end

function interactive_batch_configs()
    println("\n", repeat("═", 72))
    println("DYNAMO — INTERACTIVE FIGURE SELECTION")
    println(repeat("═", 72))
    println("Press Enter to accept a default. Multiple selections accept 1,3, 1-3, or all.")

    data_mode = prompt_one(
        "Simulation layout",
        ["custom", "standard"];
        labels = [
            "Any directory: one simulation or several simulations to compare",
            "Preconfigured VaryingMach / VaryingRes / VaryingRatio folders",
        ],
    )
    if data_mode == "custom"
        comparisons, comparison_repository = custom_comparison_selection()
    else
        comparison_repository = normalized_path(prompt_text(
            "Comparison root directory",
            COMPARISON_REPOSITORY,
        ))
        comparison_labels = [
            "$(comparison.key): $(comparison.folder) (" *
            join(comparison.simulations, ", ") * ")"
            for comparison in COMPARISONS
        ]
        comparisons = prompt_multiple(
            "Comparison groups",
            collect(COMPARISONS);
            labels = comparison_labels,
        )
    end
    inventories = print_snapshot_inventory(comparisons, comparison_repository)

    notebook_names = available_notebooks()
    notebook_labels = [
        "$(uppercasefirst(name)) — " *
        "$(length(figures_for_notebooks([name]))) registered figures"
        for name in notebook_names
    ]
    selected_notebooks = prompt_multiple(
        "Figure families",
        notebook_names;
        labels = notebook_labels,
    )
    candidate_figures = figures_for_notebooks(selected_notebooks)
    figure_labels = replace.(candidate_figures, '_' => ' ')
    selected_figures = prompt_multiple(
        "Figures to compute",
        candidate_figures;
        labels = figure_labels,
    )

    lines_of_sight = prompt_multiple(
        "Lines of sight for projected figures",
        ["x", "y", "z"],
    )
    format = prompt_one(
        "Output format",
        ["png", "pdf"];
        labels = ["PNG (recommended for quick inspection)", "PDF (vector output)"],
    )
    output_root = normalized_path(prompt_text(
        "Output root directory",
        joinpath(PROJECT_DIRECTORY, "figures"),
    ))

    available_snapshot_count = minimum(
        length(inventory.sources) for inventory in inventories
    )
    snapshot_count = prompt_positive_integer(
        "Number of snapshots per simulation (evenly spaced; Enter uses all listed)",
        available_snapshot_count,
    )

    configurations = BatchConfig[]
    for comparison in comparisons
        append_figure_jobs!(
            configurations,
            comparison,
            selected_figures,
            :even,
            snapshot_count,
            "selected",
            lines_of_sight;
            comparison_repository = comparison_repository,
            output_root = output_root,
            output_format = format,
        )
    end
    configurations, comparison_repository
end

function print_batch_plan(configurations, comparison_repository)
    println("\n", repeat("═", 72))
    println("DYNAMO — COMPARISON PLAN")
    println(repeat("═", 72))
    println("Comparison root     : ", abspath(comparison_repository))
    println("Planned jobs        : ", length(configurations))
    for (index, config) in enumerate(configurations)
        println(
            "  ",
            index,
            ". ",
            basename(normpath(config.data_repository)),
            " [",
            join(config.simulations, ", "),
            "]",
            " — ",
            config.snapshot_window,
            " ",
            config.snapshot_count,
            " snapshots — ",
            length(config.figures),
            " figures — ",
            basename(config.output_directory),
            basename(config.output_directory) == "common" ?
                " (LOS-independent, computed once)" :
                " (LOS = $(config.line_of_sight))",
        )
    end
    println(repeat("═", 72))
end

function main()
    interactive = any(argument -> argument in ("--interactive", "-i"), ARGS) ||
        lowercase(get(ENV, "DYNAMO_INTERACTIVE", "false")) in
            ("1", "true", "yes")
    configurations, comparison_repository = interactive ?
        interactive_batch_configs() : (BATCH_CONFIGS, COMPARISON_REPOSITORY)
    print_batch_plan(configurations, comparison_repository)
    lowercase(get(ENV, "DYNAMO_DRY_RUN", "false")) in ("1", "true", "yes") &&
        return

    try
        destinations = String[]
        for config in configurations
            append!(destinations, run_batch(config))
        end
        if !isempty(destinations)
            println("\nDOWNLOAD TO THE LAPTOP")
            println("Run this command from the local Git repository:")
            println("  bash download_figures.sh")
            println(
                "For future runs, compute and download in one local command:",
            )
            println("  bash run_cluster_and_download.sh")
        end
    catch error_value
        println(stderr, "\nError: ", sprint(showerror, error_value))
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
