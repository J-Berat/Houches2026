module VorticityComparison

using CairoMakie
using HDF5
using LaTeXStrings
using Serialization

export plot_vorticity_time_comparison, vorticity_time_series

const PROJECT_DIRECTORY = normpath(joinpath(@__DIR__, ".."))
const SNAPSHOT_EXTENSIONS = Set([".h5", ".hdf5"])
const PC_CM = 3.0856775814913673e18
const KM_CM = 1.0e5
const MYR_S = 3.15576e13
const CACHE_VERSION = 1

normalize_field_name(value) =
    lowercase(replace(String(value), r"[^a-zA-Z0-9]" => ""))

function snapshot_number(path)
    matched = match(r"(\d+)(?=\D*$)", splitext(basename(path))[1])
    isnothing(matched) ? typemax(Int) : parse(Int, matched.captures[1])
end

function snapshot_directory(simulation_directory)
    direct = joinpath(simulation_directory, "DataCubes")
    isdir(direct) ? direct : simulation_directory
end

"""
Return the HDF5 snapshots of one simulation.

If a `DataCubes/` directory exists it is used directly, which prevents power
spectrum HDF5 files from being mistaken for physical cubes.
"""
function discover_hdf5_snapshots(simulation_directory)
    root = snapshot_directory(abspath(expanduser(simulation_directory)))
    isdir(root) || error("Simulation directory not found: $(root)")
    files = String[]
    for (directory, subdirectories, names) in walkdir(root)
        filter!(name -> lowercase(name) ∉
            ("powerspectra", "power_spectra", "figures", "exports"),
            subdirectories)
        for name in names
            lowercase(splitext(name)[2]) in SNAPSHOT_EXTENSIONS || continue
            push!(files, joinpath(directory, name))
        end
    end
    sort!(files; by = path -> (snapshot_number(path), path))
    isempty(files) && error(
        "No HDF5 cube was found in $(root). Convert raw RAMSES outputs first " *
        "with tools/ramses_to_hdf5.py.",
    )
    files
end

function collect_dataset_paths!(paths, group, prefix = "")
    for raw_name in keys(group)
        name = String(raw_name)
        path = isempty(prefix) ? name : string(prefix, "/", name)
        object = group[name]
        try
            if object isa HDF5.Dataset
                push!(paths, path)
            elseif object isa HDF5.Group
                collect_dataset_paths!(paths, object, path)
            end
        finally
            close(object)
        end
    end
    paths
end

function dataset_path(paths, aliases; required = true)
    normalized_aliases = Set(normalize_field_name.(aliases))
    exact = findfirst(paths) do path
        normalize_field_name(basename(path)) in normalized_aliases
    end
    !isnothing(exact) && return paths[exact]
    nested = findfirst(paths) do path
        normalize_field_name(path) in normalized_aliases
    end
    !isnothing(nested) && return paths[nested]
    required && error(
        "Required HDF5 field not found. Accepted names: " *
        join(aliases, ", ") * ". Available datasets: " * join(paths, ", "),
    )
    nothing
end

function read_velocity_and_time(path; velocity_unit_kms = 1.0)
    h5open(path, "r") do handle
        paths = collect_dataset_paths!(String[], handle)
        vx_path = dataset_path(paths,
            ["vx", "velx", "velocityx", "xvelocity", "velocity/x"])
        vy_path = dataset_path(paths,
            ["vy", "vely", "velocityy", "yvelocity", "velocity/y"])
        vz_path = dataset_path(paths,
            ["vz", "velz", "velocityz", "zvelocity", "velocity/z"])
        time_path = dataset_path(paths,
            ["t", "time", "simulationtime", "myrtime"]; required = false)
        vx = Float64.(read(handle[vx_path])) .* velocity_unit_kms
        vy = Float64.(read(handle[vy_path])) .* velocity_unit_kms
        vz = Float64.(read(handle[vz_path])) .* velocity_unit_kms
        size(vx) == size(vy) == size(vz) || error(
            "Velocity components have inconsistent shapes in $(path).",
        )
        ndims(vx) == 3 || error(
            "Velocity fields must be three-dimensional in $(path).",
        )
        stored_time = if isnothing(time_path)
            NaN
        else
            raw_time = read(handle[time_path])
            Float64(raw_time isa Number ? raw_time : first(raw_time))
        end
        vx, vy, vz, stored_time
    end
end

"""
Compute volume-weighted mean and RMS vorticity without allocating three
cube-sized curl arrays. Velocities are in km s⁻¹ and the returned vorticities
are in Myr⁻¹.
"""
function vorticity_metrics(
        vx,
        vy,
        vz;
        box_length_pc = (100.0, 100.0, 100.0),
    )
    nx, ny, nz = size(vx)
    dx, dy, dz = Float64.(box_length_pc) ./ (nx, ny, nz)
    conversion = MYR_S / (PC_CM / KM_CM)
    # Julia may expose several thread pools, so threadid() can be larger than
    # nthreads() for the default pool. maxthreadid() is the safe accumulator
    # size on both Julia 1.11 and 1.12.
    accumulator_count = Threads.maxthreadid()
    thread_sum = zeros(Float64, accumulator_count)
    thread_square = zeros(Float64, accumulator_count)
    thread_count = zeros(Int, accumulator_count)

    Threads.@threads for k in 1:nz
        thread = Threads.threadid()
        km = k == 1 ? nz : k - 1
        kp = k == nz ? 1 : k + 1
        local_sum = 0.0
        local_square = 0.0
        local_count = 0
        @inbounds for j in 1:ny
            jm = j == 1 ? ny : j - 1
            jp = j == ny ? 1 : j + 1
            for i in 1:nx
                im = i == 1 ? nx : i - 1
                ip = i == nx ? 1 : i + 1
                wx = conversion * (
                    (vz[i, jp, k] - vz[i, jm, k]) / (2dy) -
                    (vy[i, j, kp] - vy[i, j, km]) / (2dz)
                )
                wy = conversion * (
                    (vx[i, j, kp] - vx[i, j, km]) / (2dz) -
                    (vz[ip, j, k] - vz[im, j, k]) / (2dx)
                )
                wz = conversion * (
                    (vy[ip, j, k] - vy[im, j, k]) / (2dx) -
                    (vx[i, jp, k] - vx[i, jm, k]) / (2dy)
                )
                magnitude = sqrt(wx^2 + wy^2 + wz^2)
                isfinite(magnitude) || continue
                local_sum += magnitude
                local_square += magnitude^2
                local_count += 1
            end
        end
        thread_sum[thread] += local_sum
        thread_square[thread] += local_square
        thread_count[thread] += local_count
    end

    count = sum(thread_count)
    count > 0 || return (mean = NaN, rms = NaN)
    (
        mean = sum(thread_sum) / count,
        rms = sqrt(sum(thread_square) / count),
    )
end

function evenly_selected(paths, maximum_snapshots)
    isnothing(maximum_snapshots) && return paths
    count = min(Int(maximum_snapshots), length(paths))
    count > 0 || error("maximum_snapshots must be positive.")
    paths[unique(round.(Int, range(1, length(paths); length = count)))]
end

function cache_path()
    joinpath(
        PROJECT_DIRECTORY,
        ".dynamo_cache",
        "vorticity_time_comparison_v$(CACHE_VERSION).bin",
    )
end

function load_cache()
    path = cache_path()
    isfile(path) || return Dict{Any,Any}()
    try
        payload = deserialize(path)
        payload isa NamedTuple || return Dict{Any,Any}()
        payload.version == CACHE_VERSION || return Dict{Any,Any}()
        payload.values isa AbstractDict || return Dict{Any,Any}()
        Dict{Any,Any}(payload.values)
    catch
        Dict{Any,Any}()
    end
end

function save_cache(values)
    path = cache_path()
    mkpath(dirname(path))
    temporary = tempname(dirname(path))
    serialize(temporary, (version = CACHE_VERSION, values))
    mv(temporary, path; force = true)
    path
end

"""
Compute the vorticity time series for one simulation directory.

All discovered snapshots are used unless `maximum_snapshots` is provided, in
which case that many snapshots are distributed evenly over the full run.
"""
function vorticity_time_series(
        simulation_directory;
        maximum_snapshots = nothing,
        box_length_pc = (100.0, 100.0, 100.0),
        velocity_unit_kms = 1.0,
        time_unit_myr = 1.0,
        label = basename(normpath(simulation_directory)),
        cache = load_cache(),
    )
    paths = evenly_selected(
        discover_hdf5_snapshots(simulation_directory),
        maximum_snapshots,
    )
    points = NamedTuple[]
    total = length(paths)
    for (index, path) in enumerate(paths)
        status = stat(path)
        key = (
            abspath(path),
            status.size,
            status.mtime,
            Tuple(Float64.(box_length_pc)),
            Float64(velocity_unit_kms),
            Float64(time_unit_myr),
        )
        print(
            "\r[",
            lpad(index, ndigits(total)),
            "/",
            total,
            "] ",
            label,
            ": ",
            basename(path),
        )
        flush(stdout)
        point = get(cache, key, nothing)
        if isnothing(point)
            vx, vy, vz, stored_time = read_velocity_and_time(
                path;
                velocity_unit_kms,
            )
            metrics = vorticity_metrics(vx, vy, vz; box_length_pc)
            fallback_time = snapshot_number(path)
            time = isfinite(stored_time) ?
                stored_time * time_unit_myr : Float64(fallback_time)
            point = (; t = time, metrics.mean, metrics.rms, path)
            cache[key] = point
        end
        push!(points, point)
    end
    println()
    sort!(points; by = point -> point.t)
    points
end

"""
Plot cooling and isothermal vorticity evolution on the same axis.

Solid lines show `⟨|ω|⟩`; dashed lines show `ω_rms`. The function saves the
figure and returns `(figure, output_path, series)`.
"""
function plot_vorticity_time_comparison(;
        isothermal_path =
            "/Xnfs/Houches2026/DynSim/isothermal_correct/turb_rms_50_N128",
        cooling_path =
            "/Xnfs/Houches2026/DynSim/cooling_freq_output/turb_rms_10_N256",
        output_path = joinpath(
            PROJECT_DIRECTORY,
            "figures",
            "cooling_isothermal_vorticity_time.png",
        ),
        maximum_snapshots = nothing,
        box_length_pc = (100.0, 100.0, 100.0),
        velocity_unit_kms = 1.0,
        time_unit_myr = 1.0,
    )
    reduction_cache = load_cache()
    specifications = [
        (
            label = "Isothermal",
            path = isothermal_path,
            color = Makie.wong_colors()[1],
        ),
        (
            label = "Cooling",
            path = cooling_path,
            color = Makie.wong_colors()[2],
        ),
    ]
    series = Dict{String,Vector{NamedTuple}}()
    for specification in specifications
        println(
            "\n",
            specification.label,
            " simulation: ",
            abspath(expanduser(specification.path)),
        )
        series[specification.label] = vorticity_time_series(
            specification.path;
            maximum_snapshots,
            box_length_pc,
            velocity_unit_kms,
            time_unit_myr,
            label = specification.label,
            cache = reduction_cache,
        )
    end
    cache_file = save_cache(reduction_cache)

    figure = Figure(size = (980, 650))
    axis = Axis(
        figure[1, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel =
            L"\langle|\omega|\rangle,\ \omega_{\mathrm{rms}}\;[\mathrm{Myr}^{-1}]",
        title = L"\mathrm{Cooling\ versus\ isothermal\ vorticity}",
        xminorticks = IntervalsBetween(5),
        yminorticks = IntervalsBetween(5),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = true,
        rightspinevisible = true,
    )
    for specification in specifications
        points = series[specification.label]
        times = getfield.(points, :t)
        means = getfield.(points, :mean)
        rms = getfield.(points, :rms)
        lines!(axis, times, means;
            color = specification.color, linewidth = 2.8)
        scatter!(axis, times, means;
            color = specification.color, markersize = 7)
        lines!(axis, times, rms;
            color = specification.color, linewidth = 2.5,
            linestyle = :dash)
        scatter!(axis, times, rms;
            color = specification.color, marker = :rect, markersize = 6)
    end
    Legend(
        figure[2, 1],
        [
            LineElement(color = specification.color, linewidth = 2.8)
            for specification in specifications
        ],
        [latexstring("\\mathrm{", specification.label, "}")
            for specification in specifications],
        L"\mathrm{Simulation}";
        orientation = :horizontal,
        framevisible = false,
        tellwidth = false,
    )
    Legend(
        figure[3, 1],
        [
            LineElement(color = :gray20, linewidth = 2.8),
            LineElement(
                color = :gray20,
                linewidth = 2.5,
                linestyle = :dash,
            ),
        ],
        [L"\langle|\omega|\rangle", L"\omega_{\mathrm{rms}}"],
        L"\mathrm{Statistic}";
        orientation = :horizontal,
        framevisible = false,
        tellwidth = false,
    )
    rowgap!(figure.layout, 5)

    destination = abspath(expanduser(output_path))
    mkpath(dirname(destination))
    save(destination, figure)
    println("\nFigure saved to: ", destination)
    println("Metric cache: ", cache_file)
    (; figure, output_path = destination, series)
end

end
