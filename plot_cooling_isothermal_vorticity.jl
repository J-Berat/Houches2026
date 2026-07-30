import Pkg

const PROJECT_DIRECTORY = @__DIR__
Pkg.activate(PROJECT_DIRECTORY)

include(joinpath(PROJECT_DIRECTORY, "src", "VorticityComparison.jl"))
using .VorticityComparison

function main(arguments = ARGS)
    isothermal_path = length(arguments) >= 1 ? arguments[1] :
        "/Xnfs/Houches2026/DynSim/isothermal_correct/turb_rms_50_N128"
    cooling_path = length(arguments) >= 2 ? arguments[2] :
        "/Xnfs/Houches2026/DynSim/cooling_freq_output/turb_rms_10_N256"
    output_path = length(arguments) >= 3 ? arguments[3] :
        joinpath(
            PROJECT_DIRECTORY,
            "figures",
            "cooling_isothermal_vorticity_time.png",
        )
    maximum_snapshots = length(arguments) >= 4 ?
        parse(Int, arguments[4]) : nothing

    plot_vorticity_time_comparison(;
        isothermal_path,
        cooling_path,
        output_path,
        maximum_snapshots,
        box_length_pc = (100.0, 100.0, 100.0),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
