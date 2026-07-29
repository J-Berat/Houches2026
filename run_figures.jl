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
        data_repository = joinpath(
            comparison_repository,
            comparison.folder,
        ),
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
    other_notebooks = filter(!=("dynamo"), SELECTED_NOTEBOOKS)
    for comparison in COMPARISONS
        comparison.key in SELECTED_COMPARISONS || continue
        "dynamo" in SELECTED_NOTEBOOKS && append_figure_jobs!(
            configurations,
            comparison,
            figures_for_notebooks(["dynamo"]),
            :first,
            20,
            "dynamo",
        )
        isempty(other_notebooks) || append_figure_jobs!(
            configurations,
            comparison,
            figures_for_notebooks(other_notebooks),
            :last,
            10,
            "observables",
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

function interactive_batch_configs()
    println("\n", repeat("═", 72))
    println("DYNAMO — INTERACTIVE FIGURE SELECTION")
    println(repeat("═", 72))
    println("Press Enter to accept a default. Multiple selections accept 1,3, 1-3, or all.")

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

    snapshot_policy = prompt_one(
        "Snapshot selection",
        ["standard", "custom"];
        labels = [
            "Scientific defaults: Dynamo first 20; observables last 10",
            "One custom snapshot window for every selected figure",
        ],
    )

    configurations = BatchConfig[]
    if snapshot_policy == "standard"
        dynamo_figures = filter(
            figure -> figure in figures_for_notebooks(["dynamo"]),
            selected_figures,
        )
        observable_figures = setdiff(selected_figures, dynamo_figures)
        for comparison in comparisons
            isempty(dynamo_figures) || append_figure_jobs!(
                configurations,
                comparison,
                dynamo_figures,
                :first,
                20,
                "dynamo",
                lines_of_sight;
                comparison_repository = comparison_repository,
                output_root = output_root,
                output_format = format,
            )
            isempty(observable_figures) || append_figure_jobs!(
                configurations,
                comparison,
                observable_figures,
                :last,
                10,
                "observables",
                lines_of_sight;
                comparison_repository = comparison_repository,
                output_root = output_root,
                output_format = format,
            )
        end
    else
        snapshot_window = Symbol(prompt_one(
            "Snapshot window",
            ["first", "last", "even"];
            labels = [
                "First snapshots",
                "Last snapshots",
                "Evenly spaced snapshots across the run",
            ],
        ))
        snapshot_count = prompt_positive_integer(
            "Number of snapshots per simulation",
            10,
        )
        for comparison in comparisons
            append_figure_jobs!(
                configurations,
                comparison,
                selected_figures,
                snapshot_window,
                snapshot_count,
                "selected",
                lines_of_sight;
                comparison_repository = comparison_repository,
                output_root = output_root,
                output_format = format,
            )
        end
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
            basename(config.data_repository),
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
