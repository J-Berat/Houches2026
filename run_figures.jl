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

# Every figure from the selected notebooks will be computed. Available groups:
# "dynamo", "dust", "starlightpol", "zeeman", "moose", and "shine".
const SELECTED_NOTEBOOKS = [
    "dynamo",
    "dust",
    "starlightpol",
    "zeeman",
    "moose",
    "shine",
]

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

function batch_config(
        comparison,
        notebooks,
        snapshot_window,
        snapshot_count,
        output_group,
    )
    BatchConfig(
        data_repository = joinpath(
            COMPARISON_REPOSITORY,
            comparison.folder,
        ),
        simulations = comparison.simulations,
        snapshot = :last,
        snapshot_window = snapshot_window,
        snapshot_count = snapshot_count,
        line_of_sight = "z",
        figures = figures_for_notebooks(notebooks),
        output_directory = joinpath(
            PROJECT_DIRECTORY,
            "figures",
            comparison.output_name,
            "$(output_group)_$(snapshot_window)$(snapshot_count)",
        ),
        output_format = "png",
    )
end

const BATCH_CONFIGS = let
    configurations = BatchConfig[]
    other_notebooks = filter(!=("dynamo"), SELECTED_NOTEBOOKS)
    for comparison in COMPARISONS
        comparison.key in SELECTED_COMPARISONS || continue
        "dynamo" in SELECTED_NOTEBOOKS && push!(
            configurations,
            batch_config(
                comparison,
                ["dynamo"],
                :first,
                20,
                "dynamo",
            ),
        )
        isempty(other_notebooks) || push!(
            configurations,
            batch_config(
                comparison,
                other_notebooks,
                :last,
                10,
                "observables",
            ),
        )
    end
    configurations
end

# =============================================================================
# EXECUTION — DO NOT EDIT BELOW THIS LINE
# =============================================================================

function print_batch_plan()
    println("\n", repeat("═", 72))
    println("DYNAMO — COMPARISON PLAN")
    println(repeat("═", 72))
    println("Comparison root     : ", abspath(COMPARISON_REPOSITORY))
    println("Planned jobs        : ", length(BATCH_CONFIGS))
    for (index, config) in enumerate(BATCH_CONFIGS)
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
            " figures",
        )
    end
    println(repeat("═", 72))
end

function main()
    print_batch_plan()
    lowercase(get(ENV, "DYNAMO_DRY_RUN", "false")) in ("1", "true", "yes") &&
        return

    try
        destinations = String[]
        for config in BATCH_CONFIGS
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
