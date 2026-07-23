import Pkg

const PROJECT_DIRECTORY = @__DIR__
Pkg.activate(PROJECT_DIRECTORY)

include(joinpath(PROJECT_DIRECTORY, "src", "DynamoAnalysis.jl"))
using .DynamoAnalysis

# =============================================================================
# CONFIGURATION TO EDIT
# =============================================================================

const DEFAULT_DATA_REPOSITORY = "/Xnfs/Houches2026/DynSim"

const CONFIG = BatchConfig(
    # Shared data root used by the notebooks on the server. The engine
    # recursively discovers simulations and their DataCubes directories.
    data_repository = get(
        ENV,
        "DYNAMO_DATA_REPOSITORY",
        DEFAULT_DATA_REPOSITORY,
    ),
    simulations = [
        "run_turb_cooling_mhd_lo_mach",
        "run_turb_cooling_mhd_mi_mach",
        "run_turb_cooling_mhd_hi_mach",
    ],
    snapshot = :last,
    line_of_sight = "z",
    figures = [
        "pdfs",
        "phase_diagram",
        "magnetic_density",
    ],
    output_directory = joinpath(PROJECT_DIRECTORY, "figures"),
    output_format = "png",
)

# =============================================================================
# EXECUTION — DO NOT EDIT BELOW THIS LINE
# =============================================================================

try
    run_batch(CONFIG)
catch error_value
    println(stderr, "\nError: ", sprint(showerror, error_value))
    exit(1)
end
