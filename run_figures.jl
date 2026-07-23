import Pkg

const PROJECT_DIRECTORY = @__DIR__
Pkg.activate(PROJECT_DIRECTORY)

include(joinpath(PROJECT_DIRECTORY, "src", "DynamoAnalysis.jl"))
using .DynamoAnalysis

# =============================================================================
# CONFIGURATION À MODIFIER
# =============================================================================

const DEFAULT_DATA_REPOSITORY = "/Xnfs/Houches2026/DynSim"

const CONFIG = BatchConfig(
    # Racine commune utilisée par les notebooks sur le serveur. Le moteur
    # découvre récursivement les simulations et leurs dossiers DataCubes.
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
# EXÉCUTION — NE RIEN MODIFIER SOUS CETTE LIGNE
# =============================================================================

try
    run_batch(CONFIG)
catch error_value
    println(stderr, "\nErreur : ", sprint(showerror, error_value))
    exit(1)
end
