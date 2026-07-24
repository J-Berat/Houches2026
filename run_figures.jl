import Pkg

const PROJECT_DIRECTORY = @__DIR__
Pkg.activate(PROJECT_DIRECTORY)

include(joinpath(PROJECT_DIRECTORY, "src", "DynamoAnalysis.jl"))
using .DynamoAnalysis

# =============================================================================
# CONFIGURATION TO EDIT
# =============================================================================

const DEFAULT_DATA_REPOSITORY = "/Xnfs/Houches2026/DynSim"

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

const SIMULATIONS = [
    "run_turb_cooling_mhd_lo_mach",
    "run_turb_cooling_mhd_mi_mach",
    "run_turb_cooling_mhd_hi_mach",
]

function batch_config(notebooks, snapshot_window, snapshot_count)
    group_name = join(notebooks, "_")
    BatchConfig(
        # Shared data root used by the notebooks on the server. The engine
        # recursively discovers simulations and their DataCubes directories.
        data_repository = get(
            ENV,
            "DYNAMO_DATA_REPOSITORY",
            DEFAULT_DATA_REPOSITORY,
        ),
        simulations = SIMULATIONS,
        snapshot = :last,
        snapshot_window = snapshot_window,
        snapshot_count = snapshot_count,
        line_of_sight = "z",
        figures = figures_for_notebooks(notebooks),
        output_directory = joinpath(
            PROJECT_DIRECTORY,
            "figures",
            "$(group_name)_$(snapshot_window)$(snapshot_count)",
        ),
        output_format = "png",
    )
end

const BATCH_CONFIGS = let
    configurations = BatchConfig[]
    "dynamo" in SELECTED_NOTEBOOKS && push!(
        configurations,
        batch_config(["dynamo"], :first, 20),
    )
    other_notebooks = filter(!=("dynamo"), SELECTED_NOTEBOOKS)
    isempty(other_notebooks) || push!(
        configurations,
        batch_config(other_notebooks, :last, 10),
    )
    configurations
end

# =============================================================================
# EXECUTION — DO NOT EDIT BELOW THIS LINE
# =============================================================================

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
