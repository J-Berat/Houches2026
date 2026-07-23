### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 34d6b5f4-9c2d-42d5-9034-543aeb8ae151
using LaTeXStrings, Printf

# ╔═╡ 8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686
begin
    using CairoMakie
    using FFTW
    using FITSIO
    using HDF5
    using PlutoUI
    using Random
    using Statistics
    using StatsBase

    K_B_CGS = 1.380649e-16          # erg K^-1
    H_PLANCK_CGS = 6.62607015e-27   # erg s
    M_H_CGS = 1.6735575e-24         # g
    C_LIGHT_CGS = 2.99792458e10     # cm s^-1
    C_LIGHT_KMS = 299_792.458        # km s^-1
    PC_CM = 3.0856775814913673e18    # cm
    MYR_S = 3.15576e13               # s
    KM_CM = 1.0e5                    # cm
    GAUSS_TO_MICROGAUSS = 1.0e6

    # Persistent caches. This cell has no reactive dependency, so Pluto runs it
    # once per session and the caches survive every widget change. RAW_CUBE_CACHE
    # holds at most one unit-free cube, exactly as stored on disk;
    # REDUCTION_CACHE holds scalar summaries, which are small enough to keep for
    # every snapshot of every run.
    const CACHE_LOCK = ReentrantLock()
    const RAW_CUBE_CACHE = Dict{Any,Any}()
    const RAW_CUBE_CACHE_BYTES = Dict{Any,Int}()
    const RAW_CUBE_CACHE_ORDER = Any[]
    const REDUCTION_CACHE = Dict{Any,Any}()
    const SNAPSHOT_SOURCE_CACHE = Dict{String,Vector{String}}()
    const DIRECTORY_DISCOVERY_CACHE = Dict{String,Vector{String}}()
    const SNAPSHOT_FINGERPRINT_CACHE = Dict{String,Any}()
    const UNREADABLE_DIRECTORY_WARNINGS = Set{String}()
    const LOCAL_HDF5_STAGE = Ref{Any}(nothing)
    const LOCAL_HDF5_STAGE_DIRECTORY = Ref{Union{Nothing,String}}(nothing)

    function positive_integer_setting(name, default)
        value = tryparse(Int, strip(get(ENV, name, string(default))))
        isnothing(value) && error("$name must be a positive integer.")
        value > 0 || error("$name must be a positive integer.")
        value
    end

    # Pluto stays deliberately conservative. The batch engine raises the entry
    # limit to the number of selected simulations, while this byte ceiling keeps
    # large production cubes from exhausting memory.
    const RAW_CUBE_CACHE_MAX_ENTRIES =
        positive_integer_setting("DYNAMO_RAW_CUBE_CACHE_ENTRIES", 1)
    const RAW_CUBE_CACHE_MAX_BYTES = positive_integer_setting(
        "DYNAMO_RAW_CUBE_CACHE_MIB",
        max(1, Int(Sys.total_memory() ÷ 4 ÷ 1024^2)),
    ) * 1024^2

    "Memoize remote file metadata for the duration of the Pluto session."
    function cached_snapshot_fingerprint(path)
        canonical = abspath(path)
        lock(CACHE_LOCK) do
            get!(SNAPSHOT_FINGERPRINT_CACHE, canonical) do
                isdir(canonical) || return (mtime(canonical), filesize(canonical))
                Tuple((basename(entry), mtime(entry), filesize(entry))
                    for entry in sort(readdir(canonical; join = true)))
            end
        end
    end

    "Return a memoized scalar summary, or `nothing` if it has not been computed."
    reduction_hit(key) = lock(() -> get(REDUCTION_CACHE, key, nothing), CACHE_LOCK)

    """
    List a directory during recursive data discovery.

    Shared simulation roots can contain unrelated directories that are not
    readable by every user. Report each such path once and skip it instead of
    aborting discovery of the accessible `DataCubes` directories.
    """
    function discovery_readdir(path; join = false)
        canonical = abspath(path)
        try
            readdir(canonical; join)
        catch error_value
            if error_value isa Base.IOError || error_value isa SystemError
                should_warn = lock(CACHE_LOCK) do
                    canonical in UNREADABLE_DIRECTORY_WARNINGS && return false
                    push!(UNREADABLE_DIRECTORY_WARNINGS, canonical)
                    true
                end
                should_warn && println(
                    stderr,
                    "Warning: skipping unreadable directory during data " *
                    "discovery: $canonical (",
                    sprint(showerror, error_value),
                    ")",
                )
                return String[]
            end
            rethrow()
        end
    end

    "Memoize a scalar summary and return it."
    store_reduction!(key, value) =
        lock(() -> (REDUCTION_CACHE[key] = value), CACHE_LOCK)

    """
    Return the cached raw cube for `key`, evaluating `build()` on a miss.

    Entries use LRU eviction and are bounded by both a count and a byte budget.
    Interactive Pluto defaults to one cube. Batch comparisons may retain several
    cubes when memory allows, preventing the same files from being reopened for
    each requested comparative figure.
    """
    function cached_raw_cube!(key, build)
        evict_oldest_unlocked!() = begin
            oldest = popfirst!(RAW_CUBE_CACHE_ORDER)
            delete!(RAW_CUBE_CACHE, oldest)
            delete!(RAW_CUBE_CACHE_BYTES, oldest)
            nothing
        end

        hit = lock(CACHE_LOCK) do
            if haskey(RAW_CUBE_CACHE, key)
                position = findfirst(isequal(key), RAW_CUBE_CACHE_ORDER)
                isnothing(position) || deleteat!(RAW_CUBE_CACHE_ORDER, position)
                push!(RAW_CUBE_CACHE_ORDER, key)
                RAW_CUBE_CACHE[key]
            else
                nothing
            end
        end
        isnothing(hit) || return hit

        # Use the largest resident cube as a conservative estimate of the next
        # allocation. Evict before reading to bound peak memory, not only the
        # final cache size.
        lock(CACHE_LOCK) do
            estimated_bytes = isempty(RAW_CUBE_CACHE_BYTES) ? 0 :
                maximum(values(RAW_CUBE_CACHE_BYTES))
            while !isempty(RAW_CUBE_CACHE_ORDER) && (
                    length(RAW_CUBE_CACHE_ORDER) >= RAW_CUBE_CACHE_MAX_ENTRIES ||
                    (estimated_bytes > 0 &&
                        sum(values(RAW_CUBE_CACHE_BYTES)) + estimated_bytes >
                            RAW_CUBE_CACHE_MAX_BYTES))
                evict_oldest_unlocked!()
            end
        end

        value = build()
        bytes = Base.summarysize(value)
        lock(CACHE_LOCK) do
            while !isempty(RAW_CUBE_CACHE_ORDER) && (
                    length(RAW_CUBE_CACHE_ORDER) >= RAW_CUBE_CACHE_MAX_ENTRIES ||
                    sum(values(RAW_CUBE_CACHE_BYTES)) + bytes >
                        RAW_CUBE_CACHE_MAX_BYTES)
                evict_oldest_unlocked!()
            end
            RAW_CUBE_CACHE[key] = value
            RAW_CUBE_CACHE_BYTES[key] = bytes
            push!(RAW_CUBE_CACHE_ORDER, key)
        end
        value
    end
    const HDF5_READ_LOCK = ReentrantLock()

    "Open at most one HDF5 file at a time and always close it before returning."
    function with_hdf5_file(operation::Function, path)
        lock(HDF5_READ_LOCK) do
            h5open(path, "r") do handle
                operation(handle)
            end
        end
    end

    """
    Whether an interactive HDF5 snapshot should first be copied to node-local
    storage. `auto` stages only `/Xnfs` paths; `true` and `false` provide an
    explicit override through `DYNAMO_LOCAL_HDF5_CACHE`.
    """
    function local_hdf5_cache_enabled(path)
        isfile(path) &&
            lowercase(splitext(path)[2]) in (".h5", ".hdf5") || return false
        setting = lowercase(strip(get(ENV, "DYNAMO_LOCAL_HDF5_CACHE", "auto")))
        setting in ("0", "false", "no", "off") && return false
        setting in ("1", "true", "yes", "on") && return true
        setting == "auto" || error(
            "DYNAMO_LOCAL_HDF5_CACHE must be auto, true, or false; found $(setting).")
        startswith(normpath(abspath(path)), "/Xnfs/")
    end

    """
    Copy `path` into a private temporary directory, retaining at most one staged
    HDF5 file. This function is called while `HDF5_READ_LOCK` is held.
    """
    function stage_hdf5_snapshot_unlocked(path)
        fingerprint = cached_snapshot_fingerprint(path)
        previous = LOCAL_HDF5_STAGE[]
        if !isnothing(previous) && previous.source == path &&
                previous.fingerprint == fingerprint && isfile(previous.staged_path)
            return previous.staged_path
        end

        if isnothing(LOCAL_HDF5_STAGE_DIRECTORY[])
            requested_parent =
                get(ENV, "DYNAMO_LOCAL_CACHE_DIRECTORY", tempdir())
            parent = requested_parent == "~" ? homedir() :
                startswith(requested_parent, "~/") ?
                    joinpath(homedir(), requested_parent[3:end]) :
                    requested_parent
            isdir(parent) || mkpath(parent)
            LOCAL_HDF5_STAGE_DIRECTORY[] =
                mktempdir(parent; prefix = "dynamo_hdf5_")
        end
        if !isnothing(previous) && isfile(previous.staged_path)
            rm(previous.staged_path; force = true)
        end
        LOCAL_HDF5_STAGE[] = nothing

        destination = joinpath(
            LOCAL_HDF5_STAGE_DIRECTORY[],
            string(hash((abspath(path), fingerprint)), "_", basename(path)),
        )
        cp(path, destination; force = true)
        LOCAL_HDF5_STAGE[] =
            (source = path, fingerprint = fingerprint, staged_path = destination)
        destination
    end

    """
    Open an analysis snapshot, optionally from the one-file node-local cache.
    The lock covers both staging and reading, so another Pluto task cannot evict
    the local file while it is in use.
    """
    function with_analysis_hdf5_file(operation::Function, path; stage_local = false)
        lock(HDF5_READ_LOCK) do
            read_path = stage_local && local_hdf5_cache_enabled(path) ?
                stage_hdf5_snapshot_unlocked(path) : path
            h5open(read_path, "r") do handle
                operation(handle)
            end
        end
    end

    """
    Render every numeric axis and colorbar tick as a genuine LaTeXString.
    """
    function latex_number(x)
        if isnan(x)
            return L"\mathrm{NaN}"
        elseif isinf(x)
            return x > 0 ? L"+\infty" : L"-\infty"
        end
        x == 0 && return L"0"
        exponent = floor(Int, log10(abs(x)))
        if exponent <= -3 || exponent >= 4
            mantissa = x / 10.0^exponent
            latexstring(@sprintf("%.3g", mantissa), "\\times 10^{", exponent, "}")
        else
            latexstring(@sprintf("%.4g", x))
        end
    end

    latex_ticklabels(values) = latex_number.(values)
    as_latex(label::LaTeXString) = label
    as_latex(label::AbstractString) = latexstring(label)

    """
    Create a Makie axis whose labels and numeric tick labels are always
    rendered by MathTeXEngine, including axes with locally customized ticks.
    """
    function latex_axis(parent; xlabel = L"", ylabel = L"",
            xtickformat = latex_ticklabels, ytickformat = latex_ticklabels,
            kwargs...)
        axis_options = (; kwargs...)
        x_is_log = get(axis_options, :xscale, identity) === log10
        y_is_log = get(axis_options, :yscale, identity) === log10
        readable_log_ticks = (
            xticksize = x_is_log ? 9 : 6,
            yticksize = y_is_log ? 9 : 6,
            xminorticksize = x_is_log ? 5 : 3,
            yminorticksize = y_is_log ? 5 : 3,
            xtickwidth = x_is_log ? 1.4 : 1.0,
            ytickwidth = y_is_log ? 1.4 : 1.0,
            xminortickwidth = x_is_log ? 1.0 : 0.8,
            yminortickwidth = y_is_log ? 1.0 : 0.8,
            xgridvisible = x_is_log,
            ygridvisible = y_is_log,
            xgridcolor = (:gray45, 0.18),
            ygridcolor = (:gray45, 0.18),
            xminorgridvisible = x_is_log,
            yminorgridvisible = y_is_log,
            xminorgridcolor = (:gray55, 0.08),
            yminorgridcolor = (:gray55, 0.08),
        )
        options = merge(readable_log_ticks, axis_options)
        Axis(parent;
            xlabel = as_latex(xlabel), ylabel = as_latex(ylabel),
            xtickformat, ytickformat, options...)
    end

    function latex_colorbar(parent, plot; label = L"",
            tickformat = latex_ticklabels, kwargs...)
        Colorbar(parent, plot;
            label = as_latex(label), tickformat, kwargs...)
    end

    """
    Render a completed Makie figure once before handing it to Pluto.

    This prevents CairoMakie's reactive text pipeline from drawing a multi-group
    legend between separate text and position updates. The original Figure is
    retained for the notebook's PNG/PDF export registry.
    """
    function stable_pluto_figure(enabled, figure)
        enabled || return nothing
        buffer = IOBuffer()
        show(buffer, MIME"image/png"(), figure)
        PlutoUI.Show(MIME"image/png"(), take!(buffer))
    end

    MHD_COLORS = [
        RGBf(0.20, 0.34, 0.84),
        RGBf(0.93, 0.47, 0.08),
        RGBf(0.04, 0.62, 0.53),
        RGBf(0.50, 0.31, 0.82),
        RGBf(0.84, 0.24, 0.48),
        RGBf(0.04, 0.57, 0.72),
    ]

    CairoMakie.activate!(type = "png")
    set_theme!(Theme(
        fontsize = 15,
        Figure = (backgroundcolor = RGBf(0.973, 0.981, 0.995),),
        Axis = (
            backgroundcolor = RGBf(1.0, 1.0, 1.0),
            xgridvisible = false, ygridvisible = false,
            topspinevisible = false, rightspinevisible = false,
            xtickformat = latex_ticklabels, ytickformat = latex_ticklabels,
        ),
        Legend = (framevisible = false,),
        Colorbar = (tickformat = latex_ticklabels, spinewidth = 0.8),
        palette = (color = MHD_COLORS,),
    ))

    notebook_style = html"""
    <style>
    :root {
        --mhd-primary: #1d4ed8;
        --mhd-secondary: #0f766e;
        --mhd-accent: #d97706;
        --mhd-violet: #7c3aed;
        --mhd-rose: #be185d;
        --mhd-success: #15803d;
        --mhd-text: #172033;
        --mhd-muted: #526079;
        --mhd-surface: #ffffff;
        --mhd-surface-soft: #f3f7fc;
        --mhd-code-surface: #eef3fb;
        --mhd-border: #cbd7e8;
        --mhd-shadow: rgba(29, 78, 216, 0.10);
        --mhd-section: var(--mhd-primary);
    }
    @media (prefers-color-scheme: dark) {
        :root {
            --mhd-primary: #79a7ff;
            --mhd-secondary: #5eead4;
            --mhd-accent: #fbbf24;
            --mhd-violet: #c4b5fd;
            --mhd-rose: #f9a8d4;
            --mhd-success: #86efac;
            --mhd-text: #edf4ff;
            --mhd-muted: #b8c5d9;
            --mhd-surface: #172033;
            --mhd-surface-soft: #1d2940;
            --mhd-code-surface: #111a2b;
            --mhd-border: #42516a;
            --mhd-shadow: rgba(0, 0, 0, 0.22);
        }
    }
    pluto-notebook {
        background: linear-gradient(160deg, var(--mhd-surface-soft) 0%, transparent 38%);
    }
    main {
        max-width: 1280px;
        padding-left: 3.5rem;
        padding-right: 3.5rem;
    }
    pluto-cell {
        --mhd-section: var(--mhd-primary);
        border-radius: 0.8rem;
    }
    pluto-cell > pluto-output {
        box-sizing: border-box;
        background: color-mix(in srgb, var(--mhd-surface) 97%, var(--mhd-section));
        border: 1px solid color-mix(in srgb, var(--mhd-border) 82%, var(--mhd-section));
        border-left: 4px solid var(--mhd-section);
        border-radius: 0.8rem;
        box-shadow: 0 5px 18px var(--mhd-shadow);
        padding: 0.7rem 1rem;
        transition: border-color 160ms ease, box-shadow 160ms ease, transform 160ms ease;
    }
    pluto-cell:hover > pluto-output {
        border-color: color-mix(in srgb, var(--mhd-section) 55%, var(--mhd-border));
        box-shadow: 0 9px 26px var(--mhd-shadow);
    }
    pluto-cell > pluto-output:has(> div:empty) {
        min-height: 0;
        background: transparent;
        border: 0;
        box-shadow: none;
        padding: 0;
    }
    pluto-cell:is(
        [id="34d6b5f4-9c2d-42d5-9034-543aeb8ae151"],
        [id="8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686"],
        [id="bdc11245-d76c-45fe-b79d-7b64861f5f53"]
    ) > pluto-output {
        min-height: 0;
        background: transparent;
        border: 0;
        box-shadow: none;
        padding: 0;
    }
    pluto-cell > pluto-input .cm-editor {
        background: var(--mhd-code-surface);
        border: 1px solid var(--mhd-border);
        border-left: 4px solid var(--mhd-violet);
        border-radius: 0.65rem;
        overflow: hidden;
    }
    pluto-cell > pluto-input .cm-gutters {
        background: color-mix(in srgb, var(--mhd-code-surface) 88%, var(--mhd-violet));
        border-right-color: var(--mhd-border);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] > pluto-output {
        background:
            radial-gradient(circle at 92% 18%, color-mix(in srgb, var(--mhd-secondary) 22%, transparent), transparent 28%),
            linear-gradient(135deg, color-mix(in srgb, var(--mhd-primary) 16%, var(--mhd-surface)), var(--mhd-surface));
        border: 1px solid color-mix(in srgb, var(--mhd-primary) 45%, var(--mhd-border));
        border-left: 7px solid var(--mhd-primary);
        border-radius: 1.05rem;
        padding: 1.3rem 1.55rem;
        box-shadow: 0 14px 38px color-mix(in srgb, var(--mhd-primary) 18%, transparent);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) ul {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 0.55rem 1rem;
        padding-left: 0;
        list-style: none;
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) li {
        background: color-mix(in srgb, var(--mhd-surface) 82%, var(--mhd-primary));
        border: 1px solid color-mix(in srgb, var(--mhd-primary) 24%, var(--mhd-border));
        border-radius: 0.55rem;
        padding: 0.65rem 0.75rem;
    }
    pluto-cell[id="3df9ad9a-b865-41f5-8d7a-34a52d0292dd"] { --mhd-section: var(--mhd-success); }
    pluto-cell[id="353cc6fb-c801-448a-a2c5-23dfd1541704"] { --mhd-section: var(--mhd-accent); }
    pluto-cell[id="e734297f-506e-45e1-8cb7-b2ae671893eb"] {
        position: sticky;
        top: 0.35rem;
        z-index: 30;
        --mhd-section: var(--mhd-accent);
    }
    pluto-cell[id="e734297f-506e-45e1-8cb7-b2ae671893eb"] > pluto-output {
        background: color-mix(in srgb, var(--mhd-surface) 94%, transparent);
        backdrop-filter: blur(12px);
        box-shadow: 0 10px 30px color-mix(in srgb, var(--mhd-accent) 18%, transparent);
        padding: 0.4rem 0.75rem;
    }
    pluto-cell[id="e734297f-506e-45e1-8cb7-b2ae671893eb"] :is(markdown, .markdown) table {
        margin: 0;
        box-shadow: none;
    }
    pluto-cell[id="e734297f-506e-45e1-8cb7-b2ae671893eb"] :is(markdown, .markdown) :is(th, td) {
        padding-top: 0.22rem;
        padding-bottom: 0.22rem;
    }

    /* Section families: setup/spatial, statistics, evolution, observations, export. */
    pluto-cell:is(
        [id="72626861-42ce-4ac0-b980-78f498f8a629"],
        [id="fa155a62-da75-4530-b5e9-215fd4f66412"],
        [id="1be45ca2-bd2f-4b77-9188-5b338d41483b"],
        [id="298bd579-bb28-48b7-8c55-ea74804b9837"],
        [id="d87379b7-3527-45a8-bc60-bec191c499af"],
        [id="a8558c31-7dcf-433e-9950-a59e9acf158b"],
        [id="24e60849-1c70-4df3-bd17-57d29949b7a6"]
    ) { --mhd-section: var(--mhd-secondary); }
    pluto-cell:is(
        [id="36aef377-3de7-435a-af83-3a90421e3159"],
        [id="dcc8f8f3-daaf-4d4b-92a3-4919ed5e36de"]
    ) { --mhd-section: var(--mhd-accent); }
    pluto-cell:is(
        [id="62440e86-b560-44ad-bb0a-43ae62e73fc3"],
        [id="d6a2f4b1-59ac-4e77-a10a-4b74c0d89231"],
        [id="6f4e2d11-2a88-41f4-93dc-01b51d86fb4f"],
        [id="67f95c39-1888-4d23-a2c2-2ee3a6cd7f0f"],
        [id="62b61ef2-8e5d-4fe9-a435-e18fb5be9461"],
        [id="ab1a7df4-ae91-47db-cb4e-cf6d42fb0143"],
        [id="478ec2f3-e057-4720-809c-17ca0a3dac21"],
        [id="14e7606a-3a13-4c8e-b860-e40dc63a6fa2"]
    ) { --mhd-section: var(--mhd-rose); }
    pluto-cell[id="3ca92f4b-1f96-4f3e-9536-ce17ae786cc6"] { --mhd-section: var(--mhd-success); }
    :is(markdown, .markdown) {
        color: var(--mhd-text);
        font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    :is(markdown, .markdown) h1 {
        color: var(--mhd-primary);
        letter-spacing: -0.025em;
        border-bottom: 3px solid var(--mhd-primary);
        padding-bottom: 0.45rem;
        margin-bottom: 1rem;
    }
    :is(markdown, .markdown) h2 {
        color: var(--mhd-section);
        background: linear-gradient(90deg, color-mix(in srgb, var(--mhd-section) 13%, transparent), transparent 78%);
        border-left: 5px solid var(--mhd-section);
        border-radius: 0.35rem;
        padding: 0.55rem 0.8rem;
        margin-top: 1rem;
    }
    :is(markdown, .markdown) h3 {
        color: var(--mhd-secondary);
        border-bottom: 1px solid var(--mhd-border);
        padding-bottom: 0.25rem;
    }
    :is(markdown, .markdown) p, :is(markdown, .markdown) li {
        color: var(--mhd-text);
    }
    :is(markdown, .markdown) strong {
        color: var(--mhd-secondary);
    }
    :is(markdown, .markdown) hr {
        border: 0;
        height: 1px;
        background: linear-gradient(90deg, var(--mhd-primary), var(--mhd-border), transparent);
        margin: 1.7rem 0 1rem;
    }
    :is(markdown, .markdown) table {
        width: 100%;
        margin: 0.45rem 0 0.65rem;
        border-collapse: separate;
        border-spacing: 0;
        overflow: hidden;
        background: var(--mhd-surface);
        border: 1px solid var(--mhd-border);
        border-radius: 0.55rem;
        box-shadow: 0 8px 22px var(--mhd-shadow);
    }
    :is(markdown, .markdown) th {
        color: var(--mhd-text);
        background: color-mix(in srgb, var(--mhd-primary) 12%, var(--mhd-surface));
        border-bottom: 2px solid var(--mhd-primary);
        padding: 0.42rem 0.55rem;
        line-height: 1.2;
    }
    :is(markdown, .markdown) td {
        border-bottom: 1px solid var(--mhd-border);
        vertical-align: middle;
        padding: 0.32rem 0.55rem;
        line-height: 1.18;
    }
    :is(markdown, .markdown) tr:last-child td {
        border-bottom: 0;
    }
    :is(markdown, .markdown) tr:nth-child(even) td {
        background: color-mix(in srgb, var(--mhd-secondary) 4%, transparent);
    }
    :is(markdown, .markdown) tr:hover td {
        background: color-mix(in srgb, var(--mhd-section) 8%, var(--mhd-surface));
    }
    :is(markdown, .markdown) blockquote {
        margin: 1rem 0;
        padding: 0.7rem 1rem;
        color: var(--mhd-text);
        background: color-mix(in srgb, var(--mhd-accent) 8%, var(--mhd-surface));
        border: 1px solid color-mix(in srgb, var(--mhd-accent) 32%, var(--mhd-border));
        border-left: 5px solid var(--mhd-accent);
        border-radius: 0.5rem;
    }
    :is(markdown, .markdown) .katex-display {
        overflow-x: auto;
        background: color-mix(in srgb, var(--mhd-section) 5%, var(--mhd-surface));
        border: 1px solid color-mix(in srgb, var(--mhd-section) 20%, var(--mhd-border));
        border-radius: 0.5rem;
        padding: 0.65rem 0.85rem;
    }
    :is(markdown, .markdown) code {
        color: var(--mhd-accent);
        background: color-mix(in srgb, var(--mhd-accent) 9%, transparent);
        border-radius: 0.25rem;
        padding: 0.08rem 0.28rem;
    }
    input[type="checkbox"], input[type="range"] {
        accent-color: var(--mhd-primary);
    }
    input[type="range"] {
        width: min(14rem, 100%);
    }

    /* Keep each plot selector close to the figure it controls. Long control
       tables become compact two-up panels; ordinary scientific tables keep
       their original layout. */
    pluto-cell:has(input, select) > pluto-output {
        padding-top: 0.5rem;
        padding-bottom: 0.5rem;
    }
    pluto-cell:has(input, select) :is(markdown, .markdown) :is(p, h3) {
        margin-top: 0.45rem;
        margin-bottom: 0.4rem;
    }
    @media (min-width: 900px) {
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(2):last-child) tbody {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(2):last-child) tbody tr {
            display: grid;
            grid-template-columns: minmax(0, 1fr) minmax(7.5rem, 0.72fr);
            align-items: center;
            border-bottom: 1px solid var(--mhd-border);
        }
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(3):last-child) tbody {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(3):last-child) tbody tr {
            display: grid;
            grid-template-columns: minmax(0, 1fr) auto auto;
            align-items: center;
            border-bottom: 1px solid var(--mhd-border);
        }
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(2):last-child) tbody tr:nth-child(odd),
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(3):last-child) tbody tr:nth-child(odd) {
            border-right: 1px solid var(--mhd-border);
        }
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(2):last-child) tbody td,
        pluto-cell:has(input, select) :is(markdown, .markdown)
        table:has(thead th:nth-child(3):last-child) tbody td {
            border-bottom: 0;
        }
    }
    select, input[type="number"], input[type="text"] {
        color: var(--mhd-text);
        background: var(--mhd-surface);
        border: 1px solid var(--mhd-border);
        border-radius: 0.45rem;
        padding: 0.28rem 0.45rem;
    }
    button {
        color: white;
        background: linear-gradient(135deg, var(--mhd-primary), color-mix(in srgb, var(--mhd-primary) 65%, var(--mhd-violet)));
        border: 0;
        border-radius: 0.5rem;
        padding: 0.42rem 0.75rem;
        box-shadow: 0 4px 12px var(--mhd-shadow);
        font-weight: 650;
    }
    button:hover {
        filter: brightness(1.08);
    }
    pluto-output img, pluto-output canvas, pluto-output svg {
        border-radius: 0.45rem;
        box-shadow: 0 9px 28px var(--mhd-shadow);
    }
    aside.plutoui-toc {
        border-left: 3px solid var(--mhd-primary);
        background: color-mix(in srgb, var(--mhd-surface) 94%, transparent);
    }

    /* Polished scientific-dashboard layer. */
    :root {
        --mhd-primary: #3157d5;
        --mhd-secondary: #0b8f80;
        --mhd-accent: #e07a16;
        --mhd-violet: #7756d8;
        --mhd-rose: #cc3f72;
        --mhd-success: #23875a;
        --mhd-text: #182238;
        --mhd-muted: #65728a;
        --mhd-surface-soft: #f1f5fc;
        --mhd-border: #dbe3f0;
        --mhd-shadow: rgba(38, 55, 95, 0.10);
    }
    pluto-notebook {
        background:
            radial-gradient(circle at 7% 4%, color-mix(in srgb, var(--mhd-primary) 11%, transparent), transparent 24rem),
            radial-gradient(circle at 92% 24%, color-mix(in srgb, var(--mhd-secondary) 9%, transparent), transparent 27rem),
            linear-gradient(180deg, #f8faff 0%, var(--mhd-surface-soft) 48%, #f8fafc 100%);
        background-attachment: fixed;
        scroll-behavior: smooth;
    }
    main {
        max-width: 1240px;
        padding-top: 1.1rem;
        padding-bottom: 5rem;
    }
    pluto-cell {
        margin-bottom: 0.42rem;
    }
    pluto-cell > pluto-output {
        background: color-mix(in srgb, var(--mhd-surface) 98%, var(--mhd-section));
        border: 1px solid color-mix(in srgb, var(--mhd-border) 88%, var(--mhd-section));
        border-left-width: 1px;
        border-radius: 1rem;
        box-shadow: 0 7px 22px rgba(38, 55, 95, 0.075);
        padding: 0.85rem 1.1rem;
    }
    pluto-cell:hover > pluto-output {
        border-color: color-mix(in srgb, var(--mhd-section) 34%, var(--mhd-border));
        box-shadow: 0 12px 34px rgba(38, 55, 95, 0.12);
        transform: translateY(-1px);
    }
    pluto-cell:is(
        [id="34d6b5f4-9c2d-42d5-9034-543aeb8ae151"],
        [id="8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686"],
        [id="bdc11245-d76c-45fe-b79d-7b64861f5f53"]
    ):hover > pluto-output,
    pluto-cell > pluto-output:has(> div:empty):hover {
        transform: none;
        box-shadow: none;
    }
    :is(markdown, .markdown) {
        font-family: "Avenir Next", Avenir, Inter, ui-sans-serif, system-ui, sans-serif;
        line-height: 1.58;
    }
    :is(markdown, .markdown) h1,
    :is(markdown, .markdown) h2,
    :is(markdown, .markdown) h3 {
        font-family: "Avenir Next", Avenir, Inter, ui-sans-serif, system-ui, sans-serif;
    }

    /* Hero: a dark, luminous cover card with compact feature tiles. */
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] > pluto-output {
        position: relative;
        overflow: hidden;
        color: #f7fbff;
        background:
            radial-gradient(circle at 88% 12%, rgba(75, 232, 210, 0.28), transparent 24rem),
            radial-gradient(circle at 8% 105%, rgba(111, 126, 255, 0.36), transparent 27rem),
            linear-gradient(128deg, #111b3b 0%, #18326a 54%, #075c67 120%);
        border: 1px solid rgba(149, 190, 255, 0.32);
        border-left: 1px solid rgba(149, 190, 255, 0.32);
        border-radius: 1.35rem;
        padding: 1.8rem 2rem 1.65rem;
        box-shadow: 0 22px 54px rgba(18, 42, 94, 0.25);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] > pluto-output::after {
        content: "MHD  •  3D DIAGNOSTICS";
        position: absolute;
        top: 1.25rem;
        right: 1.5rem;
        color: #d9fffa;
        background: rgba(255, 255, 255, 0.11);
        border: 1px solid rgba(255, 255, 255, 0.20);
        border-radius: 999px;
        padding: 0.34rem 0.68rem;
        font-size: 0.68rem;
        font-weight: 750;
        letter-spacing: 0.12em;
        backdrop-filter: blur(10px);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) h1 {
        max-width: 72%;
        color: #ffffff;
        border-bottom: 0;
        margin: 0 0 0.7rem;
        padding: 0;
        font-size: clamp(2rem, 4vw, 3.15rem);
        font-weight: 760;
        line-height: 1.06;
        letter-spacing: -0.045em;
        text-wrap: balance;
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) p,
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) li {
        color: #edf5ff;
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) ul {
        gap: 0.7rem;
        margin: 1rem 0;
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) li {
        background: rgba(255, 255, 255, 0.085);
        border: 1px solid rgba(255, 255, 255, 0.16);
        border-radius: 0.78rem;
        padding: 0.72rem 0.82rem;
        backdrop-filter: blur(8px);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) strong {
        color: #75f0dd;
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) blockquote {
        color: #f8fbff;
        background: rgba(255, 255, 255, 0.09);
        border: 1px solid rgba(255, 255, 255, 0.18);
        border-left: 4px solid #ffc56f;
        border-radius: 0.7rem;
        margin: 1rem 0;
        backdrop-filter: blur(8px);
    }
    pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) blockquote p {
        color: #f8fbff;
        margin: 0;
    }

    /* Full-bleed section ribbons inside their cards. */
    pluto-cell:has(:is(markdown, .markdown) h2) > pluto-output {
        border-top: 3px solid var(--mhd-section);
    }
    pluto-cell:has(:is(markdown, .markdown) h2) :is(markdown, .markdown) > hr:first-child {
        display: none;
    }
    :is(markdown, .markdown) h2 {
        color: var(--mhd-section);
        background:
            radial-gradient(circle at 96% 50%, color-mix(in srgb, var(--mhd-section) 15%, transparent), transparent 15rem),
            linear-gradient(90deg, color-mix(in srgb, var(--mhd-section) 13%, var(--mhd-surface)), transparent 88%);
        border: 0;
        border-bottom: 1px solid color-mix(in srgb, var(--mhd-section) 24%, var(--mhd-border));
        border-radius: 0.75rem;
        padding: 0.72rem 0.92rem;
        margin: 0 0 0.85rem;
        font-size: 1.5rem;
        font-weight: 720;
        letter-spacing: -0.025em;
    }
    :is(markdown, .markdown) h3 {
        color: color-mix(in srgb, var(--mhd-section) 82%, var(--mhd-text));
        border: 0;
        padding: 0;
        margin-top: 0.9rem;
        font-size: 1.03rem;
        font-weight: 760;
        letter-spacing: 0.015em;
    }

    /* Cleaner control panels and data tables. */
    :is(markdown, .markdown) table {
        border-color: color-mix(in srgb, var(--mhd-section) 18%, var(--mhd-border));
        border-radius: 0.8rem;
        box-shadow: 0 6px 18px rgba(38, 55, 95, 0.06);
    }
    :is(markdown, .markdown) th {
        color: #ffffff;
        background: linear-gradient(112deg, var(--mhd-section), color-mix(in srgb, var(--mhd-section) 70%, var(--mhd-violet)));
        border-bottom: 0;
        font-size: 0.82rem;
        font-weight: 740;
        letter-spacing: 0.035em;
        text-transform: uppercase;
    }
    :is(markdown, .markdown) tr:nth-child(even) td {
        background: color-mix(in srgb, var(--mhd-section) 3.5%, var(--mhd-surface));
    }
    :is(markdown, .markdown) tr:hover td {
        background: color-mix(in srgb, var(--mhd-section) 8%, var(--mhd-surface));
    }
    select, input[type="number"], input[type="text"] {
        min-height: 2rem;
        background: color-mix(in srgb, var(--mhd-surface) 96%, var(--mhd-primary));
        border-color: color-mix(in srgb, var(--mhd-primary) 22%, var(--mhd-border));
        border-radius: 0.58rem;
        box-shadow: inset 0 1px 2px rgba(35, 52, 90, 0.04);
        transition: border-color 150ms ease, box-shadow 150ms ease;
    }
    select:focus, input[type="number"]:focus, input[type="text"]:focus {
        outline: none;
        border-color: var(--mhd-primary);
        box-shadow: 0 0 0 3px color-mix(in srgb, var(--mhd-primary) 16%, transparent);
    }
    input[type="checkbox"] {
        width: 1.05rem;
        height: 1.05rem;
        vertical-align: -0.15rem;
    }
    input[type="range"] {
        accent-color: var(--mhd-section);
    }
    button {
        min-height: 2.1rem;
        border-radius: 0.62rem;
        padding: 0.46rem 0.9rem;
        box-shadow: 0 7px 17px color-mix(in srgb, var(--mhd-primary) 22%, transparent);
        transition: transform 140ms ease, filter 140ms ease, box-shadow 140ms ease;
    }
    button:hover {
        transform: translateY(-1px);
        box-shadow: 0 10px 22px color-mix(in srgb, var(--mhd-primary) 28%, transparent);
    }
    :is(markdown, .markdown) .katex-display {
        background: linear-gradient(100deg, color-mix(in srgb, var(--mhd-section) 6%, var(--mhd-surface)), var(--mhd-surface));
        border-radius: 0.72rem;
    }
    pluto-output img, pluto-output canvas, pluto-output svg {
        border-radius: 0.78rem;
        box-shadow: 0 12px 32px rgba(38, 55, 95, 0.12);
    }

    /* Glass table of contents that stays readable on long notebooks. */
    aside.plutoui-toc {
        max-height: calc(100vh - 2rem);
        overflow: auto;
        color: var(--mhd-text);
        background: color-mix(in srgb, var(--mhd-surface) 88%, transparent);
        border: 1px solid color-mix(in srgb, var(--mhd-primary) 18%, var(--mhd-border));
        border-left: 4px solid var(--mhd-primary);
        border-radius: 0.9rem;
        box-shadow: 0 14px 36px rgba(38, 55, 95, 0.14);
        backdrop-filter: blur(16px);
    }
    aside.plutoui-toc header {
        color: var(--mhd-primary);
        font-weight: 760;
    }
    aside.plutoui-toc a {
        border-radius: 0.38rem;
        transition: color 130ms ease, background 130ms ease;
    }
    aside.plutoui-toc a:hover {
        color: var(--mhd-primary);
        background: color-mix(in srgb, var(--mhd-primary) 9%, transparent);
    }
    @media (max-width: 1500px) {
        aside.plutoui-toc {
            transform: translateX(calc(100% - 3rem));
            opacity: 0.72;
            transition: transform 220ms ease, opacity 180ms ease, box-shadow 180ms ease;
        }
        aside.plutoui-toc:hover,
        aside.plutoui-toc:focus-within {
            transform: translateX(0);
            opacity: 1;
            box-shadow: 0 18px 48px rgba(38, 55, 95, 0.22);
        }
    }
    @media (prefers-color-scheme: dark) {
        :root {
            --mhd-surface: #151d31;
            --mhd-surface-soft: #0f1728;
            --mhd-code-surface: #0d1424;
            --mhd-border: #35425c;
        }
        pluto-notebook {
            background:
                radial-gradient(circle at 7% 4%, rgba(90, 119, 255, 0.16), transparent 24rem),
                radial-gradient(circle at 92% 24%, rgba(32, 196, 171, 0.12), transparent 27rem),
                linear-gradient(180deg, #0c1322 0%, #111a2b 100%);
        }
    }
    @media (max-width: 760px) {
        main { padding-left: 1rem; padding-right: 1rem; }
        :is(markdown, .markdown) table { display: block; overflow-x: auto; }
        pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) ul { grid-template-columns: 1fr; }
        pluto-cell[id="e734297f-506e-45e1-8cb7-b2ae671893eb"] { position: static; }
        pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] > pluto-output {
            padding: 1.3rem 1.15rem;
        }
        pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] > pluto-output::after {
            position: static;
            display: inline-block;
            margin-bottom: 0.85rem;
        }
        pluto-cell[id="3ef88702-bef5-4eca-a151-df97aa7ec2c4"] :is(markdown, .markdown) h1 {
            max-width: none;
        }
    }
    </style>
    """
end

# ╔═╡ 3ef88702-bef5-4eca-a151-df97aa7ec2c4
md"""
# SHINE — synthetic H I emission

21-cm radiative transfer, phase-separated H I columns, velocity moments, spectra, and RGB velocity composites.

> **Reactive mode.** Open `shine.jl` with Pluto. Selecting a repository, run, snapshot, or line of sight updates all dependent products automatically.

> **Lazy startup.** Run `run_pluto.jl`, then select `shine.jl`. Pluto starts without evaluating the expensive cells; run the result cells you need and Pluto will resolve their upstream dependencies.

All dimensional quantities are converted to the physical units shown on their axes or colorbars. Projected means are density weighted unless stated otherwise, and periodic boundaries are used for spatial operations.
"""

# ╔═╡ bdc11245-d76c-45fe-b79d-7b64861f5f53
PlutoUI.TableOfContents(title = "Notebook sections", indent = true, depth = 3, aside = true)

# ╔═╡ 7bd6f2c9-ae49-4636-a251-f526ab347125
begin
    # Keep the notebook independent from its installation directory. No data
    # repository is selected unless the user explicitly provides one here or
    # through DYNAMO_DATA_REPOSITORY (used by the non-interactive batch script).
    DEFAULT_DATA_REPOSITORY =
        strip(get(ENV, "DYNAMO_DATA_REPOSITORY", ""))
    SNAPSHOT_EXTENSIONS = (".h5", ".hdf5", ".fits", ".fit", ".fts")
    MAX_SNAPSHOTS_PER_RUN = 40
    nothing
end

# ╔═╡ 55cbdbf4-e0f2-431e-a736-09f41ab7ee75
md"""
---

## 1. Data selection and physical units

### Data repository

The selected path may be any directory containing HDF5/FITS snapshots directly,
or a parent directory containing one or more simulation folders. Directory
names and nesting are unrestricted. No cube is opened when the notebook starts:
enter a path and click **Load path** when you are ready.

| Data source | Control |
|:--|:--|
| Data path | $(@bind data_repository PlutoUI.confirm(PlutoUI.TextField(90; default = DEFAULT_DATA_REPOSITORY, placeholder = "/absolute/path/to/repository-or-DataCubes"); label = "Load path")) |
"""

# ╔═╡ 98360288-85ca-4551-bdde-c12c7a329302
begin
    snapshot_extension(path) = lowercase(splitext(path)[2])
    is_snapshot_file(path) = isfile(path) && snapshot_extension(path) in SNAPSHOT_EXTENSIONS
    is_hdf5_file(path) = isfile(path) && snapshot_extension(path) in (".h5", ".hdf5")
    is_fits_file(path) = isfile(path) && snapshot_extension(path) in (".fits", ".fit", ".fts")
    snapshot_format(path) = isdir(path) ? "FITS field directory" :
        is_fits_file(path) ? "FITS multi-extension image" : "HDF5"

    normalize_field_name(name) = lowercase(replace(strip(String(name)), r"[^A-Za-z0-9]" => ""))
    const FIELD_ALIASES = Dict(
        :rho => ["rho", "density", "massdensity"],
        :P => ["p", "pressure", "press", "thermalpressure", "gaspressure"],
        :vx => ["vx", "velx", "velocityx", "xvelocity"],
        :vy => ["vy", "vely", "velocityy", "yvelocity"],
        :vz => ["vz", "velz", "velocityz", "zvelocity"],
        :bx => ["bx", "bfieldx", "magneticfieldx"],
        :by => ["by", "bfieldy", "magneticfieldy"],
        :bz => ["bz", "bfieldz", "magneticfieldz"],
        :bx_l => ["bxl", "bxleft", "bxminus"],
        :bx_r => ["bxr", "bxright", "bxplus"],
        :by_l => ["byl", "byleft", "byminus"],
        :by_r => ["byr", "byright", "byplus"],
        :bz_l => ["bzl", "bzleft", "bzminus"],
        :bz_r => ["bzr", "bzright", "bzplus"],
        :L => ["l", "boxlength", "boxsize", "domainlength"],
        :t => ["t", "time", "snapshot_time", "simtime"],
    )

    function hdf5_dataset_paths(group, prefix = "")
        paths = String[]
        for name in keys(group)
            component = String(name)
            path = isempty(prefix) ? component : string(prefix, "/", component)
            object = group[component]
            try
                if object isa HDF5.Dataset
                    push!(paths, path)
                elseif object isa HDF5.Group
                    append!(paths, hdf5_dataset_paths(object, path))
                end
            finally
                close(object)
            end
        end
        paths
    end

    function hdf5_file_is_snapshot(path)
        isfile(path) && snapshot_extension(path) in (".h5", ".hdf5") || return false
        try
            with_hdf5_file(path) do h
                any(hdf5_dataset_paths(h)) do dataset_path
                    dataset = h[dataset_path]
                    try
                        ndims(dataset) == 3
                    finally
                        close(dataset)
                    end
                end
            end
        catch
            false
        end
    end

    function hdf5_field_path(h, field; required = true, source = "HDF5 file",
            overrides = Dict{Symbol,String}(), available_paths = nothing)
        aliases = Set(normalize_field_name.(FIELD_ALIASES[field]))
        paths = isnothing(available_paths) ? hdf5_dataset_paths(h) : available_paths
        matches = filter(paths) do path
            normalize_field_name(basename(path)) in aliases
        end
        override = get(overrides, field, "")
        if !isempty(override)
            override in paths || error(
                "The selected HDF5 mapping $(field) → $(override) does not exist in $(source). " *
                "Choose a dataset available in every snapshot.")
            return override
        elseif isempty(matches)
            required || return nothing
            available = isempty(paths) ? "(none)" : join(first(paths, min(30, length(paths))), ", ")
            error(
                "HDF5 field $(field) was not found in $(source). Accepted names: " *
                join(FIELD_ALIASES[field], ", ") *
                ". Choose it explicitly in the HDF5 field mapping table. " *
                "Available datasets: $(available)",
            )
        elseif length(matches) > 1
            error(
                "Ambiguous HDF5 field $(field) in $(source): " * join(sort(matches), ", ") *
                ". Choose the intended dataset in the HDF5 field mapping table.",
            )
        end
        only(matches)
    end

    function read_hdf5_field(h, field; required = true, source = "HDF5 file",
            overrides = Dict{Symbol,String}(), available_paths = nothing)
        path = hdf5_field_path(h, field; required, source, overrides,
            available_paths)
        isnothing(path) && return nothing
        dataset = h[path]
        try
            read(dataset)
        finally
            close(dataset)
        end
    end

    hdf5_scalar_value(value) = value isa Number ? Float64(value) : Float64(only(value))

    # Face-centred fields are averaged in their own precision, so that snapshots
    # stored as Float32 are not silently widened to Float64.
    function average_faces(lower, upper)
        half = convert(float(promote_type(eltype(lower), eltype(upper))), 0.5)
        # Reuse the lower-face buffer. The fused assignment avoids allocating a
        # third full 3-D array for every magnetic component.
        @. lower = half * (lower + upper)
        lower
    end

    """
    Cheap fingerprint that changes whenever a snapshot is rewritten on disk.

    Including it in the cache keys means the caches invalidate themselves when a
    simulation is re-run, instead of serving stale arrays.
    """
    function snapshot_fingerprint(path)
        cached_snapshot_fingerprint(path)
    end

    function centered_hdf5_magnetic_component(h, centered, left, right, source;
            overrides = Dict{Symbol,String}(), available_paths = nothing)
        direct = read_hdf5_field(h, centered; required = false, source, overrides,
            available_paths)
        isnothing(direct) || return direct
        lower = read_hdf5_field(h, left; required = false, source, overrides,
            available_paths)
        upper = read_hdf5_field(h, right; required = false, source, overrides,
            available_paths)
        (isnothing(lower) || isnothing(upper)) && error(
            "HDF5 magnetic component $(centered) in $(source) requires either " *
            "$(join(FIELD_ALIASES[centered], ", ")) or both face fields " *
            "$(join(FIELD_ALIASES[left], ", ")) and " *
            "$(join(FIELD_ALIASES[right], ", ")).",
        )
        average_faces(lower, upper)
    end

    function fits_directory_is_snapshot(path)
        isdir(path) || return false
        stems = Set(
            normalize_field_name(splitext(file)[1])
            for file in discovery_readdir(path)
            if is_fits_file(joinpath(path, file))
        )
        any(alias -> normalize_field_name(alias) in stems, FIELD_ALIASES[:rho]) &&
            all(field -> any(alias -> normalize_field_name(alias) in stems,
                FIELD_ALIASES[field]), (:vx, :vy, :vz))
    end

    function snapshot_sources(cube_directory)
        isdir(cube_directory) || return String[]
        canonical = abspath(cube_directory)
        lock(CACHE_LOCK) do
            get!(SNAPSHOT_SOURCE_CACHE, canonical) do
                fits_directory_is_snapshot(canonical) && return [canonical]
                sources = String[]
                for path in discovery_readdir(canonical; join = true)
                    if is_hdf5_file(path) || is_fits_file(path) ||
                            fits_directory_is_snapshot(path)
                        push!(sources, path)
                    end
                end
                sort(sources)
            end
        end
    end

    function expand_home(path)
        path == "~" && return homedir()
        startswith(path, "~/") ? joinpath(homedir(), path[3:end]) : path
    end

    function discover_cube_directories(path; depth = 0, max_depth = 32)
        isdir(path) || return String[]
        canonical = abspath(path)
        if depth == 0
            return lock(CACHE_LOCK) do
                get!(DIRECTORY_DISCOVERY_CACHE, canonical) do
                    discover_cube_directories_uncached(canonical; depth, max_depth)
                end
            end
        end
        discover_cube_directories_uncached(canonical; depth, max_depth)
    end

    function discover_cube_directories_uncached(path; depth = 0, max_depth = 32)
        direct_snapshots = snapshot_sources(path)
        isempty(direct_snapshots) || return [abspath(path)]
        depth >= max_depth && return String[]
        found = String[]
        for entry in sort(discovery_readdir(path))
            startswith(entry, ".") && continue
            child = joinpath(path, entry)
            isdir(child) || continue
            append!(found, discover_cube_directories_uncached(child;
                depth = depth + 1, max_depth))
        end
        unique(found)
    end

    run_directories(path) = discover_cube_directories(path)
    is_dataset_root(path) = !isempty(discover_cube_directories(path))

    function resolve_data_root(repository)
        requested = abspath(expand_home(strip(repository)))
        isdir(requested) || error("Data folder does not exist: $requested")
        is_dataset_root(requested) || error(
            "No directory containing HDF5 or FITS snapshots was found " *
            "recursively under: $requested.",
        )
        requested
    end

    function run_cube_directory(root, directory)
        isabspath(directory) && !isempty(snapshot_sources(directory)) && return directory
        candidates = [
            joinpath(root, directory, "DataCubes"),
            joinpath(root, directory),
            joinpath(directory, "DataCubes"),
            directory,
        ]
        index = findfirst(path -> !isempty(snapshot_sources(path)), candidates)
        isnothing(index) && error("No snapshots found for run folder: $directory")
        abspath(candidates[index])
    end

    function first_run_snapshot(root, directory)
        cube_directory = run_cube_directory(root, directory)
        sources = snapshot_sources(cube_directory)
        isempty(sources) && error("No HDF5 or FITS snapshots found in $cube_directory")
        first(sources)
    end

    function fits_key_value(hdu, names)
        for name in names
            value = try
                first(FITSIO.read_key(hdu, uppercase(String(name))))
            catch
                nothing
            end
            isnothing(value) || return value
        end
        nothing
    end

    function fits_hdu_name(hdu)
        value = fits_key_value(hdu, ["EXTNAME", "HDUNAME", "BTYPE", "FIELD"])
        isnothing(value) ? "" : normalize_field_name(value)
    end

    function first_fits_image_hdu(fits)
        for index in 1:length(fits)
            hdu = fits[index]
            hdu isa FITSIO.ImageHDU || continue
            data = try
                read(hdu)
            catch
                nothing
            end
            !isnothing(data) && !isempty(data) && return hdu
        end
        nothing
    end

    function fits_field_hdu(fits, field; primary_fallback = false)
        aliases = Set(normalize_field_name.(FIELD_ALIASES[field]))
        for index in 1:length(fits)
            hdu = fits[index]
            hdu isa FITSIO.ImageHDU || continue
            fits_hdu_name(hdu) in aliases && return hdu
        end
        primary_fallback ? first_fits_image_hdu(fits) : nothing
    end

    function fits_directory_field_file(directory, field)
        aliases = Set(normalize_field_name.(FIELD_ALIASES[field]))
        for file in readdir(directory; join = true)
            is_fits_file(file) || continue
            normalize_field_name(splitext(basename(file))[1]) in aliases && return file
        end
        nothing
    end

    function read_fits_field(source, field; required = true, primary_fallback = false)
        result = if isdir(source)
            file = fits_directory_field_file(source, field)
            isnothing(file) ? nothing : FITS(file) do fits
                hdu = first_fits_image_hdu(fits)
                isnothing(hdu) ? nothing : read(hdu)
            end
        else
            FITS(source) do fits
                hdu = fits_field_hdu(fits, field; primary_fallback)
                isnothing(hdu) ? nothing : read(hdu)
            end
        end
        required && isnothing(result) && error(
            "FITS field $(field) was not found in $(source). Accepted names: " *
            join(FIELD_ALIASES[field], ", "))
        result
    end

    function fits_header_scalar(source, names)
        files = isdir(source) ? filter(is_fits_file, readdir(source; join = true)) : [source]
        for file in files
            value = FITS(file) do fits
                for index in 1:length(fits)
                    candidate = fits_key_value(fits[index], names)
                    candidate isa Number && return Float64(candidate)
                end
                nothing
            end
            isnothing(value) || return value
        end
        nothing
    end

    function hdf5_scalar(root, directory, names)
        path = first_run_snapshot(root, directory)
        if is_fits_file(path) || isdir(path)
            return fits_header_scalar(path, names)
        end
        with_hdf5_file(path) do h
            candidates = vcat(names,
                ["metadata/$name" for name in names],
                ["parameters/$name" for name in names])
            for candidate in candidates
                haskey(h, candidate) || continue
                dataset = h[candidate]
                value = try
                    read(dataset)
                finally
                    close(dataset)
                end
                value isa Number && return Float64(value)
                length(value) == 1 && return Float64(first(value))
            end
            nothing
        end
    end

    function directory_parameter(directory, names)
        for name in names
            pattern = Regex("(?i)(?:^|[_-])" * name * "[_=-]?([0-9]+(?:[p.][0-9]+)?)")
            matched = match(pattern, directory)
            isnothing(matched) || return parse(Float64, replace(matched.captures[1], 'p' => '.'))
        end
        nothing
    end

    function resolution_value(root, directory, overrides = Dict{Symbol,String}())
        path = first_run_snapshot(root, directory)
        if is_fits_file(path) || isdir(path)
            dimensions = size(read_fits_field(path, :rho; primary_fallback = true))
            return all(==(first(dimensions)), dimensions) ? first(dimensions) : maximum(dimensions)
        end
        with_hdf5_file(path) do h
            density_path = hdf5_field_path(h, :rho; source = path,
                overrides)
            dataset = h[density_path]
            dimensions = try
                size(dataset)
            finally
                close(dataset)
            end
            all(==(first(dimensions)), dimensions) ? first(dimensions) : maximum(dimensions)
        end
    end

    function chi_value(root, directory)
        stored = hdf5_scalar(root, directory,
            ["chi", "Chi", "forcing_chi", "compressive_to_solenoidal_ratio"])
        isnothing(stored) || return stored
        parsed = directory_parameter(directory,
            ["chi", "forcingchi", "compressive_ratio", "compratio", "ratio"])
        isnothing(parsed) || return parsed
        name = lowercase(directory)
        occursin("solenoidal", name) && !occursin("compressive", name) && return 0.0
        occursin("compressive", name) && !occursin("solenoidal", name) && return Inf
        nothing
    end

    parameter_text(value) = isinf(value) ? "∞" : @sprintf("%.4g", value)

    run_directory_name(directory) = basename(
        lowercase(basename(directory)) == "datacubes" ? dirname(directory) : directory)

    function relative_run_name(root, directory)
        run_path = lowercase(basename(directory)) == "datacubes" ? dirname(directory) : directory
        relative = relpath(run_path, root)
        relative == "." && return basename(run_path)
        replace(relative, Base.Filesystem.path_separator => " / ")
    end

    function run_label(directory, comparison_kind, root, overrides = Dict{Symbol,String}())
        name = run_directory_name(directory)
        if comparison_kind == :resolution
            return "N = $(resolution_value(root, directory, overrides))³"
        elseif comparison_kind == :ratio
            chi = chi_value(root, directory)
            !isnothing(chi) && return "χ = $(parameter_text(chi))"
            occursin(r"_lo_(ratio|chi)$", lowercase(name)) && return "Low χ"
            occursin(r"_mi_(ratio|chi)$", lowercase(name)) && return "Intermediate χ"
            occursin(r"_hi_(ratio|chi)$", lowercase(name)) && return "High χ"
            return "χ: " * titlecase(replace(name, r"^run_turb_" => "", '_' => ' '))
        elseif comparison_kind == :mach
            occursin(r"_lo_mach$", lowercase(name)) && return "Low Mach"
            occursin(r"_mi_mach$", lowercase(name)) && return "Intermediate Mach"
            occursin(r"_hi_mach$", lowercase(name)) && return "High Mach"
        end
        comparison_kind == :folder ? relative_run_name(root, directory) :
            titlecase(replace(name, r"^run_turb_" => "", '_' => ' '))
    end

    function run_sort_key(directory, comparison_kind, root,
            overrides = Dict{Symbol,String}())
        if comparison_kind == :resolution
            return (1, Float64(resolution_value(root, directory, overrides)), directory)
        elseif comparison_kind == :ratio
            chi = chi_value(root, directory)
            !isnothing(chi) && return (1, chi, directory)
            rank = occursin(r"_lo_(ratio|chi)$", lowercase(directory)) ? 1.0 :
                occursin(r"_mi_(ratio|chi)$", lowercase(directory)) ? 2.0 :
                occursin(r"_hi_(ratio|chi)$", lowercase(directory)) ? 3.0 : 4.0
            return (2, rank, directory)
        elseif comparison_kind == :folder
            return (1, 1.0, relative_run_name(root, directory))
        end
        name = lowercase(run_directory_name(directory))
        rank = occursin(r"_lo_mach$", name) ? 1.0 :
            occursin(r"_mi_mach$", name) ? 2.0 :
            occursin(r"_hi_mach$", name) ? 3.0 : 4.0
        (1, rank, directory)
    end

    nothing
end

# ╔═╡ 6a9cfec2-2c80-4d72-94c0-cdb47aa4f046
begin
    mapping_root = resolve_data_root(data_repository)
    mapping_sources = [source
        for directory in run_directories(mapping_root)
        for source in snapshot_sources(directory)
        if is_snapshot_file(source) && !is_fits_file(source)]
    HDF5_REFERENCE_FILE = isempty(mapping_sources) ? nothing : first(mapping_sources)
    HDF5_AVAILABLE_DATASETS = isnothing(HDF5_REFERENCE_FILE) ? String[] :
        with_hdf5_file(hdf5_dataset_paths, HDF5_REFERENCE_FILE)

    function hdf5_mapping_options(field)
        aliases = Set(normalize_field_name.(FIELD_ALIASES[field]))
        matches = filter(HDF5_AVAILABLE_DATASETS) do path
            normalize_field_name(basename(path)) in aliases
        end
        automatic_label = length(matches) == 1 ? "Automatic → $(only(matches))" :
            isempty(matches) ? "No alias match — choose a dataset" :
            "Ambiguous ($(length(matches)) matches) — choose a dataset"
        options = Pair{String,String}["" => automatic_label]
        append!(options, [path => path for path in HDF5_AVAILABLE_DATASETS])
        options
    end

    mapping_widget = if isnothing(HDF5_REFERENCE_FILE)
        hdf5_rho_dataset = hdf5_P_dataset = hdf5_vx_dataset = ""
        hdf5_vy_dataset = hdf5_vz_dataset = hdf5_bx_dataset = ""
        hdf5_by_dataset = hdf5_bz_dataset = hdf5_bx_l_dataset = ""
        hdf5_bx_r_dataset = hdf5_by_l_dataset = hdf5_by_r_dataset = ""
        hdf5_bz_l_dataset = hdf5_bz_r_dataset = hdf5_L_dataset = hdf5_t_dataset = ""
        md"""
        ### Field mapping

        No HDF5 snapshot was found. FITS field names are resolved from the same alias dictionary.
        """
    else
        md"""
        ### HDF5 field mapping

        Reference file: **$(HDF5_REFERENCE_FILE)**

        Automatic mapping is used only for a **single** alias match. If the result is ambiguous or absent, explicitly choose the intended dataset. A manual choice overrides aliases.

        | Physical field | Accepted aliases | HDF5 dataset |
        |:--|:--|:--|
        | Density `ρ` | **$(join(FIELD_ALIASES[:rho], ", "))** | $(@bind hdf5_rho_dataset PlutoUI.Select(hdf5_mapping_options(:rho))) |
        | Pressure `P` | **$(join(FIELD_ALIASES[:P], ", "))** | $(@bind hdf5_P_dataset PlutoUI.Select(hdf5_mapping_options(:P))) |
        | Velocity `vₓ` | **$(join(FIELD_ALIASES[:vx], ", "))** | $(@bind hdf5_vx_dataset PlutoUI.Select(hdf5_mapping_options(:vx))) |
        | Velocity `vᵧ` | **$(join(FIELD_ALIASES[:vy], ", "))** | $(@bind hdf5_vy_dataset PlutoUI.Select(hdf5_mapping_options(:vy))) |
        | Velocity `v_z` | **$(join(FIELD_ALIASES[:vz], ", "))** | $(@bind hdf5_vz_dataset PlutoUI.Select(hdf5_mapping_options(:vz))) |
        | Magnetic `Bₓ` | **$(join(FIELD_ALIASES[:bx], ", "))** | $(@bind hdf5_bx_dataset PlutoUI.Select(hdf5_mapping_options(:bx))) |
        | Magnetic `Bᵧ` | **$(join(FIELD_ALIASES[:by], ", "))** | $(@bind hdf5_by_dataset PlutoUI.Select(hdf5_mapping_options(:by))) |
        | Magnetic `B_z` | **$(join(FIELD_ALIASES[:bz], ", "))** | $(@bind hdf5_bz_dataset PlutoUI.Select(hdf5_mapping_options(:bz))) |
        | Box length `L` | **$(join(FIELD_ALIASES[:L], ", "))** | $(@bind hdf5_L_dataset PlutoUI.Select(hdf5_mapping_options(:L))) |
        | Time `t` | **$(join(FIELD_ALIASES[:t], ", "))** | $(@bind hdf5_t_dataset PlutoUI.Select(hdf5_mapping_options(:t))) |

        #### Face-centred magnetic fields (fallback)

        These are used only when no centred `Bₓ`, `Bᵧ`, or `B_z` dataset is selected or uniquely detected.

        | Physical field | HDF5 dataset |
        |:--|:--|
        | `Bₓ` left / right | $(@bind hdf5_bx_l_dataset PlutoUI.Select(hdf5_mapping_options(:bx_l))) / $(@bind hdf5_bx_r_dataset PlutoUI.Select(hdf5_mapping_options(:bx_r))) |
        | `Bᵧ` left / right | $(@bind hdf5_by_l_dataset PlutoUI.Select(hdf5_mapping_options(:by_l))) / $(@bind hdf5_by_r_dataset PlutoUI.Select(hdf5_mapping_options(:by_r))) |
        | `B_z` left / right | $(@bind hdf5_bz_l_dataset PlutoUI.Select(hdf5_mapping_options(:bz_l))) / $(@bind hdf5_bz_r_dataset PlutoUI.Select(hdf5_mapping_options(:bz_r))) |
        """
    end

    HDF5_FIELD_OVERRIDES = Dict(
        :rho => hdf5_rho_dataset, :P => hdf5_P_dataset,
        :vx => hdf5_vx_dataset, :vy => hdf5_vy_dataset, :vz => hdf5_vz_dataset,
        :bx => hdf5_bx_dataset, :by => hdf5_by_dataset, :bz => hdf5_bz_dataset,
        :bx_l => hdf5_bx_l_dataset, :bx_r => hdf5_bx_r_dataset,
        :by_l => hdf5_by_l_dataset, :by_r => hdf5_by_r_dataset,
        :bz_l => hdf5_bz_l_dataset, :bz_r => hdf5_bz_r_dataset,
        :L => hdf5_L_dataset, :t => hdf5_t_dataset,
    )
    mapping_widget
end

# ╔═╡ 18e0d2d1-56ca-46cc-a13e-77f142612b5a
begin

    ROOT = resolve_data_root(data_repository)
    unsorted_run_dirs = run_directories(ROOT)
    comparison_kind =
        all(path -> occursin("varyingres", lowercase(path)), unsorted_run_dirs) ? :resolution :
        all(path -> occursin("varyingratio", lowercase(path)), unsorted_run_dirs) ? :ratio :
        all(path -> occursin("varyingmach", lowercase(path)), unsorted_run_dirs) ? :mach : :folder
    discovered_run_dirs = sort(unsorted_run_dirs;
        by = directory -> run_sort_key(directory, comparison_kind, ROOT,
            HDF5_FIELD_OVERRIDES))
    run_pairs = [run_label(directory, comparison_kind, ROOT,
            HDF5_FIELD_OVERRIDES) => directory
        for directory in discovered_run_dirs]
    length(unique(first.(run_pairs))) == length(run_pairs) ||
        (run_pairs = [relative_run_name(ROOT, directory) => directory
            for directory in discovered_run_dirs])
    RUN_DIRS = Dict(run_pairs)
    run_labels = first.(run_pairs)
    run_colors = Dict(label => MHD_COLORS[mod1(index, length(MHD_COLORS))]
        for (index, label) in enumerate(run_labels))
    comparison_parameter = comparison_kind == :resolution ? "grid resolution N³" :
        comparison_kind == :ratio ? "χ = Ecomp/Esol" :
        comparison_kind == :mach ? "Mach number" : "simulation folder"

    function run_label_latex_source(label)
        if startswith(label, "N = ")
            value = replace(replace(label, "N = " => ""; count = 1), "³" => "")
            return "N=" * value * "^3"
        elseif startswith(label, "χ = ")
            value = replace(replace(label, "χ = " => ""; count = 1), "∞" => "\\infty")
            return "\\chi=" * value
        elseif startswith(label, "χ: ")
            value = replace(replace(label, "χ: " => ""; count = 1), " " => "~")
            return "\\chi:\\mathrm{" * value * "}"
        elseif endswith(label, " χ")
            qualifier = replace(replace(label, r"\s*χ$" => ""), " " => "~")
            return "\\mathrm{" * qualifier * "}~\\chi"
        end
        "\\mathrm{" * replace(label, " " => "~") * "}"
    end

    latex_run_label(label) = latexstring(run_label_latex_source(label))
    # Makie legends are rebuilt reactively by Pluto.  Keep their text as one
    # atomic block: changing a LaTeXString can otherwise update its glyph blocks
    # before the matching positions and trigger a ComputePipeline length error.
    legend_run_label(label) = String(label)
    legend_rate_label(gamma; fitted = false) = string(
        fitted ? "ΓB, fit = " : "ΓB = ", round(gamma; sigdigits = fitted ? 4 : 3), " Myr⁻¹")
    latex_run_field_label(field, label) =
        latexstring(field, "\\;[", run_label_latex_source(label), "]")

    function available_snapshot_files(label)
        sources = snapshot_sources(run_cube_directory(ROOT, RUN_DIRS[label]))
        isempty(sources) && error("No HDF5 or FITS snapshots found for $label in $ROOT")
        sources
    end

    function limit_snapshot_files(sources, maximum_count)
        length(sources) <= maximum_count && return sources
        indices = unique(round.(Int,
            range(1, length(sources); length = maximum_count)))
        sources[indices]
    end

    function snapshot_time_from_filename(path)
        matched = match(r"([0-9]+(?:\.[0-9]+)?)\D*$", basename(path))
        isnothing(matched) ? NaN : parse(Float64, matched.captures[1])
    end

    function snapshot_time(path)
        if is_fits_file(path) || isdir(path)
            stored = read_fits_field(path, :t; required = false)
            !isnothing(stored) && length(stored) == 1 && return Float64(first(stored))
            header_time = fits_header_scalar(path, ["TIME", "T", "SIMTIME", "MYRTIME"])
            isnothing(header_time) && return snapshot_time_from_filename(path)
            return header_time
        end
        with_hdf5_file(path) do h
            stored = read_hdf5_field(h, :t; required = false, source = path,
                overrides = HDF5_FIELD_OVERRIDES)
            isnothing(stored) && return snapshot_time_from_filename(path)
            hdf5_scalar_value(stored)
        end
    end

    available_run_files = Dict(label => available_snapshot_files(label)
        for label in run_labels)
    run_files = Dict(label => limit_snapshot_files(
            available_run_files[label], MAX_SNAPSHOTS_PER_RUN)
        for label in run_labels)
    # Do not open every HDF5 file at startup just to populate inactive temporal
    # controls. Exact physical times are read with the selected cube or during
    # an explicitly requested temporal sweep.
    run_times = Dict(label => snapshot_time_from_filename.(run_files[label])
        for label in run_labels)
    maximum_snapshot_count = maximum(length.(values(run_files)))
    run_summary = join([
        string(label, " = ", length(run_files[label]), "/",
            length(available_run_files[label]), " snapshots")
        for label in run_labels
    ], ", ")
    nothing
end

# ╔═╡ e734297f-506e-45e1-8cb7-b2ae671893eb
md"""
| Navigation | Control |
|:--|:--|
| Run | $(@bind selected_run PlutoUI.Select(run_labels; default = run_labels[cld(length(run_labels), 2)])) |
| Simulations in comparative plots | $(@bind comparison_run_selection PlutoUI.MultiSelect(run_labels; default = [run_labels[cld(length(run_labels), 2)]])) |
| Snapshot | $(@bind selected_snapshot PlutoUI.Slider(1:length(run_files[selected_run]); default = length(run_files[selected_run]), show_value = true)) |
| Line of sight | $(@bind los_name PlutoUI.Select(["x", "y", "z"]; default = "z")) |
"""

# ╔═╡ 7174c31f-f186-48f9-b66e-29cf9a1c1fe3
md"""
| Parameter | Control |
|:--|:--|
| Adiabatic index $\gamma$ | $(@bind gamma PlutoUI.NumberField(1.0:0.01:2.0; default = 5 / 3)) |
| Mean particle mass $\mu$ [$m_{\mathrm H}$] | $(@bind mean_molecular_weight PlutoUI.NumberField(0.1:0.01:5.0; default = 1.4)) |
| Density unit [$\mathrm{g\,cm^{-3}}$ per stored unit] | $(@bind density_unit_gcm3 PlutoUI.NumberField(default = 1.0e-12)) |
| Pressure unit [$\mathrm{erg\,cm^{-3}}$ per stored unit] | $(@bind pressure_unit_ergcm3 PlutoUI.NumberField(default = 1.0)) |
| Velocity unit [$\mathrm{km\,s^{-1}}$ per stored unit] | $(@bind velocity_unit_kms PlutoUI.NumberField(default = 1.0)) |
| Magnetic-field unit [$\mathrm{G}$ per stored unit] | $(@bind magnetic_unit_G PlutoUI.NumberField(default = 1.0)) |
| Length unit [$\mathrm{pc}$ per stored unit] | $(@bind length_unit_pc PlutoUI.NumberField(default = 1.0)) |
| Time unit [$\mathrm{Myr}$ per stored unit] | $(@bind time_unit_Myr PlutoUI.NumberField(default = 1.0)) |
| PDF weighting | $(@bind pdf_weighting PlutoUI.Select(["volume", "mass"]; default = "volume")) |
| Number of bins | $(@bind nbins PlutoUI.Slider(20:5:100; default = 50, show_value = true)) |
"""

# ╔═╡ 353cc6fb-c801-448a-a2c5-23dfd1541704
begin
    requested_comparison_run_labels =
        [label for label in run_labels if label in comparison_run_selection]
    requested_open_labels = unique(vcat([selected_run], requested_comparison_run_labels))
    analysis_series_labels = requested_open_labels
    comparison_run_labels = requested_comparison_run_labels
    isempty(comparison_run_labels) && (comparison_run_labels = [selected_run])
    # Loading the selected raw cube here supplies the exact stored time and
    # primes the one-entry cache. The analysis cell below therefore reuses the
    # same arrays instead of opening this HDF5 file a second time.
    active_time_value = load_raw_cube(
        run_files[selected_run][selected_snapshot]).t * time_unit_Myr
    active_time_text = isfinite(active_time_value) ?
        string(round(active_time_value; sigdigits = 6)) : "not available"
    comparison_runs_text = join(comparison_run_labels, ", ")
    active_snapshot_format = snapshot_format(run_files[selected_run][selected_snapshot])
    Markdown.parse("""
    ### Active selection

    | Quantity | Value |
    |:--|:--|
    | Run | **$(selected_run)** |
    | Snapshot | **$(selected_snapshot)** of **$(length(run_files[selected_run]))** |
    | Input format | **$(active_snapshot_format)** |
    | Physical time | **$(active_time_text)** ``\\mathrm{Myr}`` |
    | Line of sight | **$(los_name)** |
    | HDF5 access | **Sequential; one file at a time** |
    | Time-series cubes | **Loaded only when a temporal figure is enabled** |
    | Comparative simulations | **$(comparison_runs_text)** |
    | Family comparison variable | **$(comparison_parameter)** |
    | Physical convention | ``\\gamma`` = **$(gamma)**; ``\\mu`` = **$(mean_molecular_weight)** ``m_{\\mathrm H}``; Gaussian CGS; PDF weighting = **$(pdf_weighting)** |
    """)
end

# ╔═╡ a8ef96ab-0ddd-4eb2-a216-b7d96c2a9a08
begin
    function fits_box_length(source, dimensions)
        stored = read_fits_field(source, :L; required = false)
        if !isnothing(stored)
            values = vec(Float64.(stored))
            length(values) == 1 && return fill(first(values), 3)
            length(values) >= 3 && return values[1:3]
        end
        axis_lengths = [fits_header_scalar(source, ["L$(axis)", "BOXSIZE$(axis)"])
            for axis in ("X", "Y", "Z")]
        all(value -> !isnothing(value), axis_lengths) && return Float64.(axis_lengths)
        box_length = fits_header_scalar(source, ["LBOX", "BOXSIZE", "BOXLEN"])
        !isnothing(box_length) && return fill(Float64(box_length), 3)
        pixel_sizes = [fits_header_scalar(source, ["CDELT$(axis)", "CD$(axis)_$(axis)"])
            for axis in 1:3]
        all(value -> !isnothing(value), pixel_sizes) &&
            return abs.(Float64.(pixel_sizes)) .* collect(dimensions)
        Float64.(dimensions)
    end

    function centered_fits_magnetic_component(source, centered, left, right)
        direct = read_fits_field(source, centered; required = false)
        !isnothing(direct) && return direct
        lower = read_fits_field(source, left; required = false)
        upper = read_fits_field(source, right; required = false)
        (isnothing(lower) || isnothing(upper)) && error(
            "FITS magnetic component $(centered) requires either $(centered) or both $(left)/$(right).")
        average_faces(lower, upper)
    end

    function validate_cube_shapes(source, fields)
        reference_shape = size(fields.rho)
        ndims(fields.rho) == 3 || error("Density in $(source) must be a 3-D array.")
        for field in (:P, :vx, :vy, :vz, :bx, :by, :bz)
            size(getfield(fields, field)) == reference_shape || error(
                "Field $(field) in $(source) has shape $(size(getfield(fields, field))); " *
                "expected $(reference_shape).")
        end
        nothing
    end

    """
    Read one HDF5 or FITS cube in the precision it is stored in.

    Nothing here depends on the physical-unit widgets, so the result is a pure
    function of the file and of the HDF5 field mapping, and can be cached. As
    before, the HDF5 file and every dataset handle are closed before any
    conversion happens.
    """
    function read_raw_cube(path; stage_local = false)
        if is_fits_file(path) || isdir(path)
            raw_fields = (
                rho = read_fits_field(path, :rho; primary_fallback = true),
                P = read_fits_field(path, :P),
                vx = read_fits_field(path, :vx),
                vy = read_fits_field(path, :vy),
                vz = read_fits_field(path, :vz),
                bx = centered_fits_magnetic_component(path, :bx, :bx_l, :bx_r),
                by = centered_fits_magnetic_component(path, :by, :by_l, :by_r),
                bz = centered_fits_magnetic_component(path, :bz, :bz_l, :bz_r),
            )
            validate_cube_shapes(path, raw_fields)
            return (; raw_fields...,
                L = fits_box_length(path, size(raw_fields.rho)),
                t = snapshot_time(path))
        end
        raw_fields, raw_length, raw_time =
                with_analysis_hdf5_file(path; stage_local) do h
            available_paths = hdf5_dataset_paths(h)
            raw_fields = (
                rho = read_hdf5_field(h, :rho; source = path,
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                P = read_hdf5_field(h, :P; source = path,
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                vx = read_hdf5_field(h, :vx; source = path,
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                vy = read_hdf5_field(h, :vy; source = path,
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                vz = read_hdf5_field(h, :vz; source = path,
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                bx = centered_hdf5_magnetic_component(h, :bx, :bx_l, :bx_r, path;
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                by = centered_hdf5_magnetic_component(h, :by, :by_l, :by_r, path;
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
                bz = centered_hdf5_magnetic_component(h, :bz, :bz_l, :bz_r, path;
                    overrides = HDF5_FIELD_OVERRIDES, available_paths),
            )
            raw_length = read_hdf5_field(h, :L; source = path,
                overrides = HDF5_FIELD_OVERRIDES, available_paths)
            raw_time = read_hdf5_field(h, :t; required = false, source = path,
                overrides = HDF5_FIELD_OVERRIDES, available_paths)
            (raw_fields, raw_length, raw_time)
        end
        # The HDF5 file and all dataset handles are closed before conversion or
        # scientific calculations begin.
        validate_cube_shapes(path, raw_fields)
        length_values = raw_length isa Number ? fill(Float64(raw_length), 3) :
            Float64.(vec(raw_length))
        length(length_values) == 3 || error(
            "HDF5 box length in $(path) must be a scalar or contain three values; " *
            "found size $(size(raw_length)).")
        (; raw_fields..., L = length_values,
            t = isnothing(raw_time) ? snapshot_time_from_filename(path) :
                hdf5_scalar_value(raw_time))
    end

    """
    Everything that changes the bytes read for `path`, used as a cache key.

    The HDF5 field mapping is part of the key because it selects which datasets
    the reader pulls out of the file.
    """
    raw_cube_key(path) = (String(path), snapshot_fingerprint(path),
        Tuple(sort!([string(field, "=", dataset)
            for (field, dataset) in HDF5_FIELD_OVERRIDES])))

    "Return the stored, unit-free arrays of one snapshot, reading them at most once."
    load_raw_cube(path) = cached_raw_cube!(
        raw_cube_key(path), () -> read_raw_cube(path; stage_local = true))

    "Apply a unit factor to a stored field without widening its element type."
    scale_field(A, factor) = A .* convert(float(eltype(A)), factor)

    """
    One snapshot in physical units, centred and scaled.

    The disk read is cached, so changing a unit widget only re-runs the
    multiplications below.
    """
    function scale_raw_cube(raw)
        velocity_scale = Float64(velocity_unit_kms)
        magnetic_scale = Float64(magnetic_unit_G)
        (
            rho = scale_field(raw.rho, density_unit_gcm3),
            P = scale_field(raw.P, pressure_unit_ergcm3),
            vx = scale_field(raw.vx, velocity_scale),
            vy = scale_field(raw.vy, velocity_scale),
            vz = scale_field(raw.vz, velocity_scale),
            bx = scale_field(raw.bx, magnetic_scale),
            by = scale_field(raw.by, magnetic_scale),
            bz = scale_field(raw.bz, magnetic_scale),
            L = Float64(length_unit_pc) .* raw.L,
            t = Float64(time_unit_Myr) * raw.t,
        )
    end

    load_cube(path) = scale_raw_cube(load_raw_cube(path))

    """
    How many snapshots a sweep may process at once.

    Each in-flight snapshot holds a raw cube plus its scaled copy, so the worker
    count is bounded by free memory as well as by the available threads.
    """
    function sweep_concurrency(cube_bytes)
        cube_bytes > 0 || return 1
        # Sys.free_memory() is not usable here: on macOS it counts only genuinely
        # free pages and is almost always near zero, which would quietly reduce
        # the sweep to one worker. Budget a fraction of total memory instead.
        # A worker holds the raw cube, its scaled copy and the derived fields,
        # which together come to roughly three times the raw size.
        budget = Sys.total_memory() ÷ 4
        affordable = Int(budget ÷ (3 * cube_bytes))
        clamp(affordable, 1, Threads.nthreads())
    end

    "Apply `work` to every index, running at most `workers` of them at a time."
    function parallel_foreach(work, indices, workers)
        workers <= 1 && return foreach(work, indices)
        queue = Channel{eltype(indices)}(length(indices))
        foreach(index -> put!(queue, index), indices)
        close(queue)
        @sync for _ in 1:workers
            Threads.@spawn for index in queue
                work(index)
            end
        end
        nothing
    end

    "Every widget value that changes the physical content of a loaded cube."
    unit_signature() = (Float64(density_unit_gcm3), Float64(pressure_unit_ergcm3),
        Float64(velocity_unit_kms), Float64(magnetic_unit_G),
        Float64(length_unit_pc), Float64(time_unit_Myr))

    "Cache key identifying a snapshot together with the units it is loaded in."
    cube_signature(path) = (raw_cube_key(path), unit_signature())

    number_density(rho) = rho ./ (Float64(mean_molecular_weight) * M_H_CGS)

    # Reductions over a cube run in one pass and accumulate in Float64, so a
    # Float32 snapshot is never materialized as a Float64 copy. Only the
    # quantile-based helpers need the selected values collected into a vector,
    # and they now allocate that vector once instead of twice.
    function collect_finite(A, keep)
        values = Vector{Float64}(undef, length(A))
        selected = 0
        @inbounds for element in A
            value = Float64(element)
            if isfinite(value) && keep(value)
                selected += 1
                values[selected] = value
            end
        end
        resize!(values, selected)
    end

    finite_values(A) = collect_finite(A, _ -> true)
    finite_positive_values(A) = collect_finite(A, >(0))

    function finite_mean(A; default = NaN)
        total = 0.0
        counted = 0
        @inbounds for element in A
            value = Float64(element)
            if isfinite(value)
                total += value
                counted += 1
            end
        end
        counted == 0 ? default : total / counted
    end

    finite_quantile(A, q; default = NaN) = begin
        values = finite_values(A)
        isempty(values) ? default : quantile(values, q)
    end

    function finite_extrema(A; default = (NaN, NaN))
        low, high = Inf, -Inf
        @inbounds for element in A
            value = Float64(element)
            if isfinite(value)
                value < low && (low = value)
                value > high && (high = value)
            end
        end
        low > high ? default : (low, high)
    end

    function finite_sum(A)
        total = 0.0
        @inbounds for element in A
            value = Float64(element)
            isfinite(value) && (total += value)
        end
        total
    end

    function finite_maximum(A; default = NaN)
        high = -Inf
        @inbounds for element in A
            value = Float64(element)
            isfinite(value) && value > high && (high = value)
        end
        isfinite(high) ? high : default
    end

    safe_log10(x) = begin
        value = Float64(x)
        isfinite(value) && value > 0 ? log10(value) : NaN
    end

    function finite_sum_dims(A, d)
        dropdims(sum(ifelse.(isfinite.(A), Float64.(A), 0.0); dims = d); dims = d)
    end

    function finite_mean_dims(A, d)
        valid = isfinite.(A)
        numerator = dropdims(sum(ifelse.(valid, Float64.(A), 0.0); dims = d); dims = d)
        denominator = dropdims(sum(valid; dims = d); dims = d)
        map((value, count) -> count > 0 ? value / count : NaN, numerator, denominator)
    end

    function finite_maximum_dims(A, d)
        maximum_values = dropdims(maximum(ifelse.(isfinite.(A), Float64.(A), -Inf);
            dims = d); dims = d)
        map(value -> isfinite(value) ? value : NaN, maximum_values)
    end

    function symlog10(A)
        finite_nonzero = filter(x -> isfinite(x) && x > 0, abs.(vec(Float64.(A))))
        linthresh = isempty(finite_nonzero) ? 1.0 : quantile(finite_nonzero, 0.10)
        sign.(A) .* log10.(1 .+ abs.(A) ./ max(linthresh, eps()))
    end

    function robust_colorrange(A, percentile; diverging = false)
        values = finite_values(A)
        isempty(values) && return (-1.0, 1.0)
        if diverging
            limit = quantile(abs.(values), percentile / 100)
            limit = limit > 0 ? limit : 1.0
            return (-limit, limit)
        end
        tail = (1 - percentile / 100) / 2
        lo, hi = quantile(values, (tail, 1 - tail))
        if lo == hi
            delta = max(abs(lo), 1.0) * 1e-6
            lo, hi = lo - delta, hi + delta
        end
        (lo, hi)
    end

    function turbulent_velocity(c)
        T = float(eltype(c.rho))
        rho, vx, vy, vz = c.rho, c.vx, c.vy, c.vz
        wsum, sx, sy, sz = 0.0, 0.0, 0.0, 0.0
        @inbounds for index in eachindex(rho)
            weight = Float64(rho[index])
            (isfinite(weight) && weight > 0) || continue
            ux, uy, uz = Float64(vx[index]), Float64(vy[index]), Float64(vz[index])
            (isfinite(ux) && isfinite(uy) && isfinite(uz)) || continue
            wsum += weight
            sx += weight * ux
            sy += weight * uy
            sz += weight * uz
        end
        wsum > 0 || return (
            vbar = (NaN, NaN, NaN),
            dvx = fill(T(NaN), size(rho)), dvy = fill(T(NaN), size(rho)),
            dvz = fill(T(NaN), size(rho)),
            dv2 = fill(T(NaN), size(rho)),
        )
        vbar = (sx / wsum, sy / wsum, sz / wsum)
        dvx, dvy, dvz = vx .- T(vbar[1]), vy .- T(vbar[2]), vz .- T(vbar[3])
        dv2 = dvx .^ 2 .+ dvy .^ 2 .+ dvz .^ 2
        (; vbar, dvx, dvy, dvz, dv2)
    end

    function magnetic_fields(c)
        B2 = c.bx .^ 2 .+ c.by .^ 2 .+ c.bz .^ 2
        (; B2, B = sqrt.(B2))
    end

    function weighted_project(A, rho, d)
        valid = isfinite.(A) .& isfinite.(rho) .& (rho .> 0)
        weights = ifelse.(valid, Float64.(rho), 0.0)
        numerator = dropdims(sum(weights .* ifelse.(valid, Float64.(A), 0.0);
            dims = d); dims = d)
        denominator = dropdims(sum(weights; dims = d); dims = d)
        map((value, weight) -> weight > 0 ? value / weight : NaN,
            numerator, denominator)
    end

    function periodic_derivative(A, d, dx)
        forward = ntuple(i -> i == d ? -1 : 0, ndims(A))
        backward = ntuple(i -> i == d ? 1 : 0, ndims(A))
        (circshift(A, forward) .- circshift(A, backward)) ./ (2dx)
    end

    function vorticity(c)
        dx = c.L ./ size(c.rho)
        kms_per_pc_to_Myr_inv = MYR_S / (PC_CM / KM_CM)
        wx = kms_per_pc_to_Myr_inv .* (periodic_derivative(c.vz, 2, dx[2]) .- periodic_derivative(c.vy, 3, dx[3]))
        wy = kms_per_pc_to_Myr_inv .* (periodic_derivative(c.vx, 3, dx[3]) .- periodic_derivative(c.vz, 1, dx[1]))
        wz = kms_per_pc_to_Myr_inv .* (periodic_derivative(c.vy, 1, dx[1]) .- periodic_derivative(c.vx, 2, dx[2]))
        (; wx, wy, wz, magnitude = sqrt.(wx .^ 2 .+ wy .^ 2 .+ wz .^ 2))
    end

    struct DecadeTicks end

    function decade_tick_values(vmin, vmax)
        lower, upper = minmax(Float64(vmin), Float64(vmax))
        lower = max(lower, floatmin(Float64))
        upper = max(upper, lower * (1 + eps(Float64)))
        first_exponent = floor(Int, log10(lower))
        last_exponent = ceil(Int, log10(upper))
        exponent_step = max(1, cld(last_exponent - first_exponent, 12))
        exponents = collect(first_exponent:exponent_step:last_exponent)
        last(exponents) == last_exponent || push!(exponents, last_exponent)
        10.0 .^ exponents
    end

    function latex_decade_number(value)
        exponent = round(Int, log10(value))
        if -3 <= exponent <= 6
            exponent >= 0 && return latexstring(@sprintf("%.0f", value))
            return latexstring(@sprintf("%.*f", -exponent, value))
        end
        latexstring("10^{", exponent, "}")
    end

    function enclosing_decade_limits(values)
        valid = Float64.(filter(value -> isfinite(value) && value > 0, vec(values)))
        isempty(valid) && return nothing
        lower, upper = extrema(valid)
        decade_lower = 10.0^floor(Int, log10(lower))
        decade_upper = 10.0^ceil(Int, log10(upper))
        decade_lower == decade_upper && (decade_upper *= 10.0)
        # Keep the bounding major ticks slightly inside the plotted frame.
        # Makie can then construct the nine logarithmic minor ticks even when
        # the data span contains only one strictly interior decade.
        padding = 10.0^0.035
        decade_lower / padding, decade_upper * padding
    end

    Makie.get_tickvalues(::DecadeTicks, ::typeof(log10), vmin, vmax) =
        decade_tick_values(vmin, vmax)

    function Makie.get_ticks(::DecadeTicks, ::typeof(log10), formatter, vmin, vmax)
        values = decade_tick_values(vmin, vmax)
        values, latex_decade_number.(values)
    end

    DECADE_TICKS = DecadeTicks()

    function polarization_column_figure(column, fraction, ylabel, color)
        N = vec(Float64.(column))
        p_percent = 100 .* vec(Float64.(fraction))
        valid = isfinite.(N) .& isfinite.(p_percent) .& (N .> 0) .&
            (p_percent .>= 0)
        N, p_percent = N[valid], p_percent[valid]
        fig = Figure(size = (620, 430))
        if isempty(N)
            Label(fig[1, 1], L"\mathrm{No\ valid\ polarization\ samples.}";
                fontsize = 20)
            return fig
        end
        axis = latex_axis(fig[1, 1],
            xlabel = L"N_{\mathrm H}\;[\mathrm{cm}^{-2}]",
            ylabel = ylabel, xscale = log10,
            xticks = DECADE_TICKS,
            xminorticks = IntervalsBetween(9), xminorticksvisible = true)
        sample_step = max(1, cld(length(N), 6000))
        sample = 1:sample_step:length(N)
        scatter!(axis, N[sample], p_percent[sample];
            color = (color, 0.22), markersize = 5)
        log_limits = extrema(log10.(N))
        if log_limits[1] < log_limits[2]
            edges = 10 .^ range(log_limits...; length = 18)
            centers, lower, medians, upper = Float64[], Float64[], Float64[], Float64[]
            for bin in 1:length(edges)-1
                members = (N .>= edges[bin]) .& (N .< edges[bin + 1])
                count(members) >= 2 || continue
                values = p_percent[members]
                push!(centers, sqrt(edges[bin] * edges[bin + 1]))
                push!(lower, quantile(values, 0.16))
                push!(medians, median(values))
                push!(upper, quantile(values, 0.84))
            end
            if !isempty(centers)
                band!(axis, centers, lower, upper; color = (color, 0.18))
                lines!(axis, centers, medians; color, linewidth = 3)
                scatter!(axis, centers, medians; color, markersize = 7)
            end
        end
        fig
    end

    function periodic_bilinear(map, xmap, ymap)
        nx, ny = size(map)
        output = similar(xmap, Float64)
        for index in eachindex(output)
            x = mod(Float64(xmap[index]) - 1, nx) + 1
            y = mod(Float64(ymap[index]) - 1, ny) + 1
            x0, y0 = floor(Int, x), floor(Int, y)
            x1, y1 = mod1(x0 + 1, nx), mod1(y0 + 1, ny)
            fx, fy = x - x0, y - y0
            output[index] = (1 - fx) * (1 - fy) * map[x0, y0] +
                fx * (1 - fy) * map[x1, y0] +
                (1 - fx) * fy * map[x0, y1] + fx * fy * map[x1, y1]
        end
        output
    end

    function lic_texture(vx, vy; niter = 1, len = 12, normalize_vectors = true,
            amplitude_weight = true, amplitude_floor = 0.15, seed = 42)
        nx, ny = size(vx)
        magnitude = sqrt.(Float64.(vx) .^ 2 .+ Float64.(vy) .^ 2)
        safe_magnitude = max.(magnitude, eps(Float64))
        if normalize_vectors
            ux, uy = vx ./ safe_magnitude, vy ./ safe_magnitude
        else
            scale = max(maximum(safe_magnitude), eps(Float64))
            ux, uy = vx ./ scale, vy ./ scale
        end
        x0 = repeat(reshape(collect(1.0:nx), nx, 1), 1, ny)
        y0 = repeat(reshape(collect(1.0:ny), 1, ny), nx, 1)
        texture = [mod(sin((i + seed) * 12.9898 + (j + seed) * 78.233) *
            43758.5453, 1.0) for i in 1:nx, j in 1:ny]
        integration_length = max(Int(len), 1)
        for _ in 1:max(Int(niter), 1)
            source = copy(texture)
            forward_x, forward_y = copy(x0), copy(y0)
            backward_x, backward_y = copy(x0), copy(y0)
            accumulated = zeros(Float64, nx, ny)
            for _ in 0:integration_length
                forward_dx = periodic_bilinear(ux, forward_x, forward_y)
                forward_dy = periodic_bilinear(uy, forward_x, forward_y)
                backward_dx = periodic_bilinear(ux, backward_x, backward_y)
                backward_dy = periodic_bilinear(uy, backward_x, backward_y)
                forward_x .+= 0.25 .* forward_dx
                forward_y .+= 0.25 .* forward_dy
                backward_x .-= 0.25 .* backward_dx
                backward_y .-= 0.25 .* backward_dy
                accumulated .+= periodic_bilinear(source, forward_x, forward_y) .+
                    periodic_bilinear(source, backward_x, backward_y)
            end
            texture .= accumulated ./ (2 * (integration_length + 1))
        end
        if amplitude_weight
            floor_value = clamp(Float64(amplitude_floor), 0.0, 1.0)
            texture .*= floor_value .+ (1 - floor_value) .* magnitude ./
                max(maximum(magnitude), eps(Float64))
        end
        texture[magnitude .<= eps(Float64)] .= 0.0
        finite_texture = filter(isfinite, vec(texture))
        isempty(finite_texture) && return zeros(Float64, nx, ny)
        low, high = quantile(finite_texture, (0.01, 0.99))
        high > low || return zeros(Float64, nx, ny)
        clamp.((texture .- low) ./ (high - low), 0.0, 1.0)
    end
end

# ╔═╡ 32110739-b60e-4592-856a-dd74f7a37401
begin
    selected_path = run_files[selected_run][selected_snapshot]
    cube = load_cube(selected_path)
    mag = magnetic_fields(cube)
    turb = turbulent_velocity(cube)
    los_dim = Dict("x" => 1, "y" => 2, "z" => 3)[los_name]
    sky_dims = filter(!=(los_dim), (1, 2, 3))
    axis_names = ("x", "y", "z")
    sky_labels = axis_names[collect(sky_dims)]
    comparison_snapshot_indices = Dict(label =>
        min(Int(selected_snapshot), length(run_files[label]))
        for label in comparison_run_labels)
    function comparison_cube(label)
        label == selected_run && comparison_snapshot_indices[label] == selected_snapshot &&
            return cube
        load_cube(run_files[label][comparison_snapshot_indices[label]])
    end
end

# ╔═╡ 94a0a0dc-baf6-4e62-a51e-dc6124d98fd4
begin
    Bcomponents = (cube.bx, cube.by, cube.bz)
    dx_los_pc = cube.L[los_dim] / size(cube.rho, los_dim)
    dx_los_cm = dx_los_pc * PC_CM
    column_density = finite_sum_dims(cube.rho, los_dim) .* dx_los_cm ./
        (Float64(mean_molecular_weight) * M_H_CGS)
    T = Float64(mean_molecular_weight) * M_H_CGS .* cube.P ./ (K_B_CGS .* cube.rho)
    Tmean = weighted_project(T, cube.rho, los_dim)
    Blos = GAUSS_TO_MICROGAUSS .* weighted_project(Bcomponents[los_dim], cube.rho, los_dim)
    Bsky1 = GAUSS_TO_MICROGAUSS .* weighted_project(Bcomponents[sky_dims[1]], cube.rho, los_dim)
    Bsky2 = GAUSS_TO_MICROGAUSS .* weighted_project(Bcomponents[sky_dims[2]], cube.rho, los_dim)
    Bsky = sqrt.(Bsky1 .^ 2 .+ Bsky2 .^ 2)
    Bmean_projected = GAUSS_TO_MICROGAUSS .* weighted_project(mag.B, cube.rho, los_dim)
    velocity_projected = weighted_project(sqrt.(turb.dv2), cube.rho, los_dim)
    magnetic_energy_density = mag.B2 ./ (8pi)
    kinetic_energy_density = 0.5 .* cube.rho .* turb.dv2 .* KM_CM^2
    thermal_energy_density = gamma > 1 + sqrt(eps(Float64)) ? cube.P ./ (gamma - 1) : cube.P
    project_energy_density(A) = finite_sum_dims(A, los_dim) .* dx_los_cm
    magnetic_energy_map = project_energy_density(magnetic_energy_density)
    kinetic_energy_map = project_energy_density(kinetic_energy_density)
    thermal_energy_map = project_energy_density(thermal_energy_density)
    omega = vorticity(cube)
    omega_map = finite_mean_dims(omega.magnitude, los_dim)
    enstrophy = 0.5 .* omega.magnitude .^ 2
    enstrophy_map = finite_mean_dims(enstrophy, los_dim)
    sky_coordinates = (
        range(0, cube.L[sky_dims[1]]; length = size(cube.rho, sky_dims[1])),
        range(0, cube.L[sky_dims[2]]; length = size(cube.rho, sky_dims[2])),
    )
end

# ╔═╡ c12d1f54-40b8-4865-9562-8dcb519f924a
md"""
---

## 2. Projected spatial diagnostics

### Heatmap selection

**Display projected maps:** $(@bind display_projected_maps PlutoUI.CheckBox(default = true))

| Field | Display | Logarithmic scale |
|:--|:--:|:--:|
| Column density | $(@bind show_column PlutoUI.CheckBox(default = true)) | $(@bind log_column PlutoUI.CheckBox(default = true)) |
| Mean temperature | $(@bind show_temperature PlutoUI.CheckBox(default = true)) | $(@bind log_temperature PlutoUI.CheckBox(default = false)) |
| Mean $B_{\mathrm{LOS}}$ | $(@bind show_blos PlutoUI.CheckBox(default = true)) | $(@bind log_blos PlutoUI.CheckBox(default = false)) *(symmetric log)* |
| Mean $B_{\mathrm{sky}}$ | $(@bind show_bsky PlutoUI.CheckBox(default = true)) | $(@bind log_bsky PlutoUI.CheckBox(default = false)) |
| Mean magnetic-field strength | $(@bind show_Bmean PlutoUI.CheckBox(default = false)) | $(@bind log_Bmean PlutoUI.CheckBox(default = false)) |
| Mean turbulent velocity | $(@bind show_velocity PlutoUI.CheckBox(default = false)) | $(@bind log_velocity PlutoUI.CheckBox(default = false)) |
| Mean vorticity | $(@bind show_vorticity PlutoUI.CheckBox(default = false)) | $(@bind log_vorticity PlutoUI.CheckBox(default = true)) |
| Integrated magnetic energy | $(@bind show_magnetic_energy PlutoUI.CheckBox(default = false)) | $(@bind log_magnetic_energy PlutoUI.CheckBox(default = true)) |
| Integrated kinetic energy | $(@bind show_kinetic_energy PlutoUI.CheckBox(default = false)) | $(@bind log_kinetic_energy PlutoUI.CheckBox(default = true)) |
| Integrated thermal energy | $(@bind show_thermal_energy PlutoUI.CheckBox(default = false)) | $(@bind log_thermal_energy PlutoUI.CheckBox(default = true)) |

| Magnetic-field overlay and rendering | Control |
|:--|:--|
| Projected-field arrows over column density | $(@bind show_projected_B PlutoUI.CheckBox(default = true)) |
| Arrow stride | $(@bind arrow_stride PlutoUI.Slider(2:1:16; default = 4, show_value = true)) |
| LIC texture over column density | $(@bind show_projected_B_lic PlutoUI.CheckBox(default = false)) |
| LIC integration length [pixels] | $(@bind lic_length PlutoUI.Slider(2:2:40; default = 12, show_value = true)) |
| LIC iterations | $(@bind lic_iterations PlutoUI.Slider(1:1:4; default = 1, show_value = true)) |
| Normalize field vectors before LIC | $(@bind lic_normalize_vectors PlutoUI.CheckBox(default = true)) |
| Weight LIC by $B_{\rm sky}$ amplitude | $(@bind lic_amplitude_weight PlutoUI.CheckBox(default = true)) |
| LIC amplitude floor | $(@bind lic_amplitude_floor PlutoUI.Slider(0.0:0.05:1.0; default = 0.15, show_value = true)) |
| LIC overlay opacity | $(@bind lic_opacity PlutoUI.Slider(0.05:0.05:1.0; default = 0.45, show_value = true)) |
| LIC texture seed | $(@bind lic_seed PlutoUI.NumberField(0:1:100000; default = 42)) |
| Heatmap contrast percentile | $(@bind color_percentile PlutoUI.Slider(90.0:0.5:100.0; default = 99.0, show_value = true)) |

"""

# ╔═╡ 62440e86-b560-44ad-bb0a-43ae62e73fc3
md"""
---

## 12. Shared observational beam

Optional elliptical Gaussian beam applied to $I$, $Q$, and $U$ before computing polarization.

| Gaussian PSF setting | Control |
|:--|:--|
| Apply Gaussian beam | $(@bind apply_observational_beam PlutoUI.CheckBox(default = false)) |
| Beam-width unit | $(@bind observational_beam_unit PlutoUI.Select(["Sky pixels" => "pixel", "Parsecs" => "pc"]; default = "pixel")) |
| Major-axis FWHM | $(@bind observational_beam_fwhm_major PlutoUI.NumberField(0.1:0.1:1000.0; default = 3.0)) |
| Minor-axis FWHM | $(@bind observational_beam_fwhm_minor PlutoUI.NumberField(0.1:0.1:1000.0; default = 3.0)) |
| Position angle [$^\circ$] | $(@bind observational_beam_pa_deg PlutoUI.Slider(0.0:1.0:180.0; default = 0.0, show_value = true)) |
| Display observable structure functions | $(@bind display_observational_structure_functions PlutoUI.CheckBox(default = true)) |
| Observable structure-function order $p$ | $(@bind observational_structure_order PlutoUI.Slider(1:4; default = 2, show_value = true)) |
| Number of observable separation samples | $(@bind observational_structure_samples PlutoUI.Slider(4:2:20; default = 10, show_value = true)) |
"""

# ╔═╡ 47b786d6-c7b5-44f4-946a-b8c485ad6380
begin
    function observational_beam_width_pixels(c, plane_dims)
        major = max(Float64(observational_beam_fwhm_major), eps(Float64))
        minor = max(Float64(observational_beam_fwhm_minor), eps(Float64))
        if observational_beam_unit == "pc"
            dx = c.L[plane_dims[1]] / size(c.rho, plane_dims[1])
            dy = c.L[plane_dims[2]] / size(c.rho, plane_dims[2])
            reference_pixel = sqrt(dx * dy)
            major /= reference_pixel
            minor /= reference_pixel
        end
        max(major, minor), min(major, minor)
    end

    function gaussian_beam_transfer(map_size, fwhm_major_pix, fwhm_minor_pix, pa_deg)
        nx, ny = map_size
        sigma_major = Float64(fwhm_major_pix) / (2sqrt(2log(2)))
        sigma_minor = Float64(fwhm_minor_pix) / (2sqrt(2log(2)))
        angle = deg2rad(Float64(pa_deg))
        cosine, sine = cos(angle), sin(angle)
        frequency_x = reshape(Float64.(FFTW.fftfreq(nx, 1.0)), nx, 1)
        frequency_y = reshape(Float64.(FFTW.fftfreq(ny, 1.0)), 1, ny)
        frequency_major = cosine .* frequency_x .+ sine .* frequency_y
        frequency_minor = -sine .* frequency_x .+ cosine .* frequency_y
        exp.(-2pi^2 .* (sigma_major^2 .* frequency_major .^ 2 .+
            sigma_minor^2 .* frequency_minor .^ 2))
    end

    function apply_gaussian_beam_2d(image, fwhm_major_pix, fwhm_minor_pix, pa_deg)
        transfer = gaussian_beam_transfer(size(image), fwhm_major_pix,
            fwhm_minor_pix, pa_deg)
        valid = isfinite.(image)
        all(valid) && begin
            filtered = ifft(fft(image) .* transfer)
            return eltype(image) <: Real ? real.(filtered) : filtered
        end
        filled = ifelse.(valid, image, zero(eltype(image)))
        filtered = ifft(fft(filled) .* transfer)
        normalization = real.(ifft(fft(Float64.(valid)) .* transfer))
        output = eltype(image) <: Real ? real.(filtered) : filtered
        map((value, weight) -> weight > sqrt(eps(Float64)) ? value / weight :
            convert(eltype(output), NaN), output, normalization)
    end

    function apply_observational_beam_2d(image, c, plane_dims)
        apply_observational_beam || return copy(image)
        major, minor = observational_beam_width_pixels(c, plane_dims)
        apply_gaussian_beam_2d(image, major, minor, observational_beam_pa_deg)
    end

    function apply_observational_beam_cube(cube_xyν, c, plane_dims)
        apply_observational_beam || return copy(cube_xyν)
        major, minor = observational_beam_width_pixels(c, plane_dims)
        output = similar(cube_xyν, promote_type(eltype(cube_xyν), Float64))
        @views for channel in axes(cube_xyν, 3)
            output[:, :, channel] .= apply_gaussian_beam_2d(
                cube_xyν[:, :, channel], major, minor, observational_beam_pa_deg)
        end
        output
    end

    "Axis-averaged periodic scalar structure function of a two-dimensional map."
    function scalar_structure_function_2d(field, lags, order; period = nothing)
        ndims(field) == 2 || error("Projected structure functions require a 2-D map.")
        data = Float64.(field)
        shifted = similar(data)
        values = zeros(Float64, length(lags))
        for (lag_index, lag) in pairs(lags)
            moment_sum = 0.0
            moment_count = 0
            for dimension in 1:2
                shift = ntuple(d -> d == dimension ? lag : 0, 2)
                circshift!(shifted, data, shift)
                for index in eachindex(data)
                    increment = shifted[index] - data[index]
                    if !isnothing(period)
                        increment = mod(increment + period / 2, period) - period / 2
                    end
                    moment = abs(increment)^order
                    if isfinite(moment)
                        moment_sum += moment
                        moment_count += 1
                    end
                end
            end
            values[lag_index] = moment_count > 0 ? moment_sum / moment_count : NaN
        end
        values
    end

    "Plot projected structure functions for a collection of observable maps."
    function observational_structure_figure(specs, c, plane_dims, order, samples;
            heading = "Projected observable structure functions")
        maximum_lag = max(1, minimum(size(first(specs).data)) ÷ 2)
        lags = unique(round.(Int, exp.(range(log(1.0), log(Float64(maximum_lag));
            length = Int(samples)))))
        pixel_scale_pc = minimum(c.L[dimension] / size(c.rho, dimension)
            for dimension in plane_dims)
        separations_pc = lags .* pixel_scale_pc
        ncols = min(2, length(specs))
        nrows = cld(length(specs), ncols)
        figure = Figure(size = (540ncols, 390nrows + 55))
        Label(figure[0, 1:ncols], heading; fontsize = 22, font = :bold)
        for (spec_index, spec) in enumerate(specs)
            row, column = cld(spec_index, ncols), mod1(spec_index, ncols)
            axis = latex_axis(figure[row, column],
                xlabel = L"\ell\;[\mathrm{pc}]",
                ylabel = latexstring("S_{", order, "}(\\ell)"),
                title = as_latex(spec.label), xscale = log10, yscale = log10,
                xticks = DECADE_TICKS, yticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                xminorticksvisible = true, yminorticksvisible = true)
            values = scalar_structure_function_2d(spec.data, lags, order;
                period = spec.period)
            valid = isfinite.(values) .& (values .> 0)
            if any(valid)
                lines!(axis, separations_pc[valid], values[valid];
                    color = spec.color, linewidth = 2.5)
                scatter!(axis, separations_pc[valid], values[valid];
                    color = spec.color, markersize = 6)
            else
                text!(axis, 0.5, 0.5; text = "constant or invalid map",
                    space = :relative, align = (:center, :center))
            end
        end
        figure
    end

    function moose_instrument_transfer(map_size, largest_scale_pix, smallest_scale_pix)
        nx, ny = map_size
        largest = max(Float64(largest_scale_pix), eps(Float64))
        smallest = max(Float64(smallest_scale_pix), 2.0)
        frequency_x = reshape(Float64.(FFTW.fftfreq(nx, 1.0)), nx, 1)
        frequency_y = reshape(Float64.(FFTW.fftfreq(ny, 1.0)), 1, ny)
        frequency2 = frequency_x .^ 2 .+ frequency_y .^ 2
        low_frequency = 1 / largest
        high_frequency = min(1 / smallest, 0.5)
        Float64.((frequency2 .>= low_frequency^2) .&
            (frequency2 .<= high_frequency^2))
    end

    function apply_moose_interferometer_2d(image, transfer)
        size(image) == size(transfer) || error("MOOSE transfer-mask shape mismatch.")
        finite_image = ifelse.(isfinite.(image), image, zero(eltype(image)))
        filtered = ifft(fft(finite_image) .* transfer)
        eltype(image) <: Real ? real.(filtered) : filtered
    end

    function apply_moose_interferometer_cube(cube_xyν, transfer)
        output = similar(cube_xyν, promote_type(eltype(cube_xyν), Float64))
        @views for channel in axes(cube_xyν, 3)
            output[:, :, channel] .= apply_moose_interferometer_2d(
                cube_xyν[:, :, channel], transfer)
        end
        output
    end

    function add_moose_qu_noise!(Q, U, snr, rng)
        signal_to_noise = Float64(snr)
        isfinite(signal_to_noise) && signal_to_noise > 0 ||
            error("MOOSE Q/U signal-to-noise ratio must be finite and positive.")
        if ndims(Q) == 2
            polarized_rms = sqrt(finite_mean(abs2.(Q); default = 0.0) +
                finite_mean(abs2.(U); default = 0.0))
            sigma = polarized_rms / signal_to_noise
            sigma > 0 && (Q .+= sigma .* randn(rng, size(Q));
                U .+= sigma .* randn(rng, size(U)))
        else
            @views for channel in axes(Q, 3)
                Qchannel, Uchannel = Q[:, :, channel], U[:, :, channel]
                polarized_rms = sqrt(finite_mean(abs2.(Qchannel); default = 0.0) +
                    finite_mean(abs2.(Uchannel); default = 0.0))
                sigma = polarized_rms / signal_to_noise
                sigma > 0 && (Qchannel .+= sigma .* randn(rng, size(Qchannel));
                    Uchannel .+= sigma .* randn(rng, size(Uchannel)))
            end
        end
        Q, U
    end
end

# ╔═╡ 478ec2f3-e057-4720-809c-17ca0a3dac21
md"""
---

## 18. SHINE H I post-processing

Synthetic $\mathrm{H\,I}$ 21-cm transfer with CNM, LNM, and WNM components.

**Display SHINE H I maps:** $(@bind display_shine PlutoUI.CheckBox(default = true))  
**Display the selected H I spectrum:** $(@bind display_shine_spectrum PlutoUI.CheckBox(default = true))  
**Display the H I velocity RGB composite:** $(@bind display_shine_rgb PlutoUI.CheckBox(default = true))

| SHINE setting | Control |
|:--|:--|
| Neutral $\mathrm{H\,I}$ fraction | $(@bind shine_neutral_fraction PlutoUI.NumberField(0.0:0.01:1.0; default = 1.0)) |
| Mass per H nucleon [$m_H$] | $(@bind shine_mu_H PlutoUI.NumberField(1.0:0.01:2.0; default = 1.4)) |
| Thermal broadening mass $\mu$ [$m_H$] | $(@bind shine_thermal_mu PlutoUI.NumberField(0.1:0.1:5.0; default = 1.0)) |
| Fixed line width [$\mathrm{km\,s^{-1}}$; $0$ uses thermal broadening] | $(@bind shine_fixed_width_kms PlutoUI.NumberField(0.0:0.1:30.0; default = 0.0)) |
| Minimum velocity [$\mathrm{km\,s^{-1}}$] | $(@bind shine_velocity_min PlutoUI.NumberField(-200.0:0.5:0.0; default = -30.0)) |
| Maximum velocity [$\mathrm{km\,s^{-1}}$] | $(@bind shine_velocity_max PlutoUI.NumberField(0.0:0.5:200.0; default = 30.0)) |
| Channel width [$\mathrm{km\,s^{-1}}$] | $(@bind shine_velocity_step PlutoUI.NumberField(0.25:0.25:10.0; default = 1.0)) |
| CNM/LNM boundary [$\mathrm{K}$] | $(@bind shine_TCNM PlutoUI.NumberField(10.0:10.0:2000.0; default = 200.0)) |
| LNM/WNM boundary [$\mathrm{K}$] | $(@bind shine_TWNM PlutoUI.NumberField(100.0:50.0:20000.0; default = 2000.0)) |
| FFT CNM threshold [$\mathrm{(km\,s^{-1})^{-1}}$] | $(@bind shine_fft_klim PlutoUI.NumberField(0.0:0.01:2.0; default = 0.12)) |
| First sky-axis pixel | $(@bind shine_sky_i PlutoUI.Slider(1:size(cube.rho, sky_dims[1]); default = cld(size(cube.rho, sky_dims[1]), 2), show_value = true)) |
| Second sky-axis pixel | $(@bind shine_sky_j PlutoUI.Slider(1:size(cube.rho, sky_dims[2]); default = cld(size(cube.rho, sky_dims[2]), 2), show_value = true)) |
| Logarithmic column-density maps | $(@bind log_shine_column PlutoUI.CheckBox(default = true)) |
| Total $N_{\mathrm{HI}}$ map | $(@bind show_shine_NHI PlutoUI.CheckBox(default = true)) |
| CNM column-density map | $(@bind show_shine_NCNM PlutoUI.CheckBox(default = false)) |
| LNM column-density map | $(@bind show_shine_NLNM PlutoUI.CheckBox(default = false)) |
| WNM column-density map | $(@bind show_shine_NWNM PlutoUI.CheckBox(default = false)) |
| Peak $T_B$ map | $(@bind show_shine_peakTb PlutoUI.CheckBox(default = true)) |
| Moment-0 map | $(@bind show_shine_mom0 PlutoUI.CheckBox(default = true)) |
| Moment-1 map | $(@bind show_shine_mom1 PlutoUI.CheckBox(default = true)) |
| Moment-2 map | $(@bind show_shine_mom2 PlutoUI.CheckBox(default = true)) |
| Peak optical-depth map | $(@bind show_shine_tau PlutoUI.CheckBox(default = false)) |
| FFT CNM tracer map | $(@bind show_shine_fftcnm PlutoUI.CheckBox(default = false)) |
| Optically thick $T_B(v)$ | $(@bind show_shine_Tb_spectrum PlutoUI.CheckBox(default = true)) |
| Optical-depth spectrum $\tau(v)$ | $(@bind show_shine_tau_spectrum PlutoUI.CheckBox(default = true)) |
| Blue/green velocity boundary [$\mathrm{km\,s^{-1}}$] | $(@bind shine_rgb_blue_green PlutoUI.NumberField(-200.0:0.5:200.0; default = 4.0)) |
| Green/red velocity boundary [$\mathrm{km\,s^{-1}}$] | $(@bind shine_rgb_green_red PlutoUI.NumberField(-200.0:0.5:200.0; default = 10.0)) |
| RGB stretch percentile | $(@bind shine_rgb_percentile PlutoUI.NumberField(90.0:0.5:100.0; default = 99.5)) |
| RGB asinh softening | $(@bind shine_rgb_softening PlutoUI.NumberField(0.01:0.01:1.0; default = 0.10)) |
| Normalize each velocity band independently | $(@bind shine_rgb_independent PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ bcc05889-02bb-47cf-b672-139e8efe4137
begin
    function shine_hi_spectrum(n, velocity, temperature, channels, dz_cm, thermal_mu, fixed_width)
        Tb = zeros(Float64, length(channels))
        tau_front = zeros(Float64, length(channels))
        for cell in eachindex(n)
            (n[cell] <= 0 || temperature[cell] <= 0) && continue
            sigma = fixed_width > 0 ? fixed_width :
                sqrt(K_B_CGS * temperature[cell] / (M_H_CGS * thermal_mu)) / KM_CM
            sigma = max(sigma, eps(Float64))
            normalization = 1 / (sqrt(2pi) * sigma)
            for channel in eachindex(channels)
                argument = (channels[channel] - velocity[cell]) / sigma
                profile = exp(-0.5 * argument^2) * normalization
                tau_cell = profile * n[cell] * dz_cm / (1.823e18 * temperature[cell])
                Tb[channel] += temperature[cell] * (-expm1(-tau_cell)) * exp(-tau_front[channel])
                tau_front[channel] += tau_cell
            end
        end
        Tb, tau_front
    end

    shine_permutation = (sky_dims[1], sky_dims[2], los_dim)
    shine_nHI_native = Float64(shine_neutral_fraction) .* cube.rho ./
        (max(Float64(shine_mu_H), eps(Float64)) * M_H_CGS)
    shine_nHI = permutedims(shine_nHI_native, shine_permutation)
    shine_temperature = permutedims(max.(T, eps(Float64)), shine_permutation)
    shine_velocity = permutedims((cube.vx, cube.vy, cube.vz)[los_dim], shine_permutation)
    shine_TCNM_value = min(Float64(shine_TCNM), Float64(shine_TWNM))
    shine_TWNM_value = max(Float64(shine_TCNM), Float64(shine_TWNM))
    shine_nCNM = shine_nHI .* (shine_temperature .< shine_TCNM_value)
    shine_nLNM = shine_nHI .* ((shine_temperature .>= shine_TCNM_value) .&
        (shine_temperature .< shine_TWNM_value))
    shine_nWNM = shine_nHI .* (shine_temperature .>= shine_TWNM_value)
    shine_NHI = finite_sum_dims(shine_nHI, 3) .* dx_los_cm
    shine_NCNM = finite_sum_dims(shine_nCNM, 3) .* dx_los_cm
    shine_NLNM = finite_sum_dims(shine_nLNM, 3) .* dx_los_cm
    shine_NWNM = finite_sum_dims(shine_nWNM, 3) .* dx_los_cm
    shine_NHI = apply_observational_beam_2d(shine_NHI, cube, sky_dims)
    shine_NCNM = apply_observational_beam_2d(shine_NCNM, cube, sky_dims)
    shine_NLNM = apply_observational_beam_2d(shine_NLNM, cube, sky_dims)
    shine_NWNM = apply_observational_beam_2d(shine_NWNM, cube, sky_dims)

    shine_vlo = min(Float64(shine_velocity_min), Float64(shine_velocity_max))
    shine_vhi = max(Float64(shine_velocity_min), Float64(shine_velocity_max))
    shine_dv = max(Float64(shine_velocity_step), 0.25)
    shine_velocity_axis = collect(shine_vlo:shine_dv:shine_vhi)
    length(shine_velocity_axis) >= 2 || (shine_velocity_axis = [shine_vlo, shine_vlo + shine_dv])
    shine_Tb = zeros(Float64, size(shine_nHI, 1), size(shine_nHI, 2), length(shine_velocity_axis))
    shine_tau = similar(shine_Tb)
    shine_mu_value = max(Float64(shine_thermal_mu), eps(Float64))
    shine_fixed_width_value = max(Float64(shine_fixed_width_kms), 0.0)
    Threads.@threads for i in axes(shine_nHI, 1)
        for j in axes(shine_nHI, 2)
            Tb_line, tau_line = shine_hi_spectrum(
                @view(shine_nHI[i, j, :]), @view(shine_velocity[i, j, :]),
                @view(shine_temperature[i, j, :]), shine_velocity_axis, dx_los_cm,
                shine_mu_value, shine_fixed_width_value)
            shine_Tb[i, j, :] .= Tb_line
            shine_tau[i, j, :] .= tau_line
        end
    end
    shine_Tb = apply_observational_beam_cube(shine_Tb, cube, sky_dims)
    shine_tau = apply_observational_beam_cube(shine_tau, cube, sky_dims)

    shine_velocity_3d = reshape(shine_velocity_axis, 1, 1, :)
    shine_positive_weight = max.(shine_Tb, 0.0) .* shine_dv
    shine_mom0 = finite_sum_dims(shine_positive_weight, 3)
    shine_weight_sum = max.(shine_mom0, eps(Float64))
    shine_mom1 = finite_sum_dims(shine_positive_weight .* shine_velocity_3d, 3) ./ shine_weight_sum
    shine_mom2 = sqrt.(max.(finite_sum_dims(shine_positive_weight .*
        (shine_velocity_3d .- reshape(shine_mom1, size(shine_mom1)..., 1)) .^ 2, 3) ./
        shine_weight_sum, 0.0))
    shine_peakTb = finite_maximum_dims(shine_Tb, 3)
    shine_peak_tau = finite_maximum_dims(shine_tau, 3)
    shine_fft_frequency = fftfreq(length(shine_velocity_axis), shine_dv)
    shine_fft_indices = findall(>(Float64(shine_fft_klim)), shine_fft_frequency)
    shine_fftcnm = zeros(Float64, size(shine_NHI))
    if !isempty(shine_fft_indices)
        Threads.@threads for i in axes(shine_Tb, 1)
            for j in axes(shine_Tb, 2)
                amplitudes = abs.(fft(@view shine_Tb[i, j, :]))
                shine_fftcnm[i, j] = amplitudes[1] > 0 ? maximum(amplitudes[shine_fft_indices]) / amplitudes[1] : 0.0
            end
        end
    end

    shine_specs = NamedTuple[]
    function add_shine_column!(enabled, data, symbol)
        enabled || return
        shown = log_shine_column ? safe_log10.(data) : data
        label = log_shine_column ? latexstring("\\log_{10}(", symbol,
            "/[\\mathrm{cm}^{-2}])") : latexstring(symbol, "\\;[\\mathrm{cm}^{-2}]")
        push!(shine_specs, (data = shown, label, colormap = :magma, diverging = false))
    end
    add_shine_column!(show_shine_NHI, shine_NHI, "N_{\\mathrm{HI}}")
    add_shine_column!(show_shine_NCNM, shine_NCNM, "N_{\\mathrm{CNM}}")
    add_shine_column!(show_shine_NLNM, shine_NLNM, "N_{\\mathrm{LNM}}")
    add_shine_column!(show_shine_NWNM, shine_NWNM, "N_{\\mathrm{WNM}}")
    show_shine_peakTb && push!(shine_specs, (data = shine_peakTb,
        label = L"\max_v T_B\;[\mathrm{K}]", colormap = :thermal, diverging = false))
    show_shine_mom0 && push!(shine_specs, (data = shine_mom0,
        label = L"M_0\;[\mathrm{K\,km\,s}^{-1}]", colormap = :viridis, diverging = false))
    show_shine_mom1 && push!(shine_specs, (data = shine_mom1,
        label = L"M_1\;[\mathrm{km\,s}^{-1}]", colormap = :balance, diverging = true))
    show_shine_mom2 && push!(shine_specs, (data = shine_mom2,
        label = L"M_2\;[\mathrm{km\,s}^{-1}]", colormap = :plasma, diverging = false))
    show_shine_tau && push!(shine_specs, (data = shine_peak_tau,
        label = L"\max_v\tau_{21}", colormap = :inferno, diverging = false))
    show_shine_fftcnm && push!(shine_specs, (data = shine_fftcnm,
        label = L"f_{\mathrm{CNM}}^{\mathrm{FFT}}", colormap = :viridis, diverging = false))

    if isempty(shine_specs)
        fig_shine = Figure(size = (900, 180))
        Label(fig_shine[1, 1], L"\mathrm{Select\ at\ least\ one\ SHINE\ H\,I\ map.}", fontsize = 20)
    else
        shine_ncols = length(shine_specs) == 1 ? 1 : 2
        shine_nrows = cld(length(shine_specs), shine_ncols)
        fig_shine = Figure(size = (550shine_ncols, 420shine_nrows))
        for (index, spec) in enumerate(shine_specs)
            row, col = cld(index, shine_ncols), mod1(index, shine_ncols)
            panel = fig_shine[row, col] = GridLayout()
            ax = latex_axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap,
                colorrange = robust_colorrange(spec.data, color_percentile; diverging = spec.diverging))
            scatter!(ax, [sky_coordinates[1][Int(shine_sky_i)]],
                [sky_coordinates[2][Int(shine_sky_j)]];
                marker = :cross, markersize = 18, strokewidth = 3, color = :white)
            latex_colorbar(panel[1, 2], hm; label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_shine ? fig_shine : nothing
end

# ╔═╡ 5ad762e1-105f-4cf7-9cf1-e0bb8c6f1bf5
begin
    "Integrate positive H I brightness over a closed velocity interval."
    function shine_velocity_band(Tb, velocities, lower, upper, dv;
            include_lower = true, include_upper = true)
        indices = findall(velocities) do velocity
            lower_ok = include_lower ? lower <= velocity : lower < velocity
            upper_ok = include_upper ? velocity <= upper : velocity < upper
            lower_ok && upper_ok
        end
        isempty(indices) && return zeros(Float64, size(Tb, 1), size(Tb, 2))
        dropdims(sum(max.(@view(Tb[:, :, indices]), 0.0); dims = 3); dims = 3) .* dv
    end

    "Apply a robust asinh display stretch to one velocity-integrated map."
    function shine_rgb_stretch(data, scale, softening)
        safe_scale = max(scale, eps(Float64))
        safe_softening = max(softening, eps(Float64))
        stretched = asinh.(max.(data, 0.0) ./ (safe_softening * safe_scale)) ./
            asinh(1 / safe_softening)
        clamp.(stretched, 0.0, 1.0)
    end

    shine_rgb_boundary_1 = clamp(min(Float64(shine_rgb_blue_green),
        Float64(shine_rgb_green_red)), shine_vlo, shine_vhi)
    shine_rgb_boundary_2 = clamp(max(Float64(shine_rgb_blue_green),
        Float64(shine_rgb_green_red)), shine_vlo, shine_vhi)

    shine_rgb_blue_map = shine_velocity_band(shine_Tb, shine_velocity_axis,
        shine_vlo, shine_rgb_boundary_1, shine_dv; include_upper = false)
    shine_rgb_green_map = shine_velocity_band(shine_Tb, shine_velocity_axis,
        shine_rgb_boundary_1, shine_rgb_boundary_2, shine_dv)
    shine_rgb_red_map = shine_velocity_band(shine_Tb, shine_velocity_axis,
        shine_rgb_boundary_2, shine_vhi, shine_dv; include_lower = false)

    shine_rgb_quantile = clamp(Float64(shine_rgb_percentile) / 100, 0.0, 1.0)
    shine_rgb_scales = [quantile(filter(isfinite, vec(channel)), shine_rgb_quantile)
        for channel in (shine_rgb_red_map, shine_rgb_green_map, shine_rgb_blue_map)]
    if !shine_rgb_independent
        shared_scale = maximum(shine_rgb_scales)
        shine_rgb_scales .= shared_scale
    end
    shine_rgb_softening_value = max(Float64(shine_rgb_softening), 0.01)
    shine_rgb_red = shine_rgb_stretch(shine_rgb_red_map, shine_rgb_scales[1],
        shine_rgb_softening_value)
    shine_rgb_green = shine_rgb_stretch(shine_rgb_green_map, shine_rgb_scales[2],
        shine_rgb_softening_value)
    shine_rgb_blue = shine_rgb_stretch(shine_rgb_blue_map, shine_rgb_scales[3],
        shine_rgb_softening_value)
    shine_rgb_image = RGBf.(shine_rgb_red, shine_rgb_green, shine_rgb_blue)

    fig_shine_rgb = Figure(size = (820, 670))
    shine_rgb_axis = latex_axis(fig_shine_rgb[1, 1],
        xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
        ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"),
        title = L"\mathrm{H\,I\ velocity\ composite}")
    image!(shine_rgb_axis,
        (first(sky_coordinates[1]), last(sky_coordinates[1])),
        (first(sky_coordinates[2]), last(sky_coordinates[2])),
        shine_rgb_image)
    rowsize!(fig_shine_rgb.layout, 1, Aspect(1, 1))

    shine_rgb_labels = [
        latexstring("R:\\ ", @sprintf("%.4g", shine_rgb_boundary_2),
            "<v_{\\mathrm{LOS}}\\leq ", @sprintf("%.4g", shine_vhi),
            "\\;\\mathrm{km\\,s}^{-1}"),
        latexstring("G:\\ ", @sprintf("%.4g", shine_rgb_boundary_1),
            "\\leq v_{\\mathrm{LOS}}\\leq ", @sprintf("%.4g", shine_rgb_boundary_2),
            "\\;\\mathrm{km\\,s}^{-1}"),
        latexstring("B:\\ ", @sprintf("%.4g", shine_vlo),
            "\\leq v_{\\mathrm{LOS}}<", @sprintf("%.4g", shine_rgb_boundary_1),
            "\\;\\mathrm{km\\,s}^{-1}"),
    ]
    shine_rgb_label_colors = [RGBf(0.86, 0.16, 0.20), RGBf(0.05, 0.62, 0.34), RGBf(0.12, 0.35, 0.92)]
    shine_rgb_legend = GridLayout(fig_shine_rgb[2, 1])
    for index in eachindex(shine_rgb_labels)
        Label(shine_rgb_legend[1, index], shine_rgb_labels[index];
            color = shine_rgb_label_colors[index], fontsize = 16, tellwidth = true)
    end
    display_shine_rgb ? fig_shine_rgb : nothing
end

# ╔═╡ 27e51ba5-4592-4766-9dde-0de383a889a0
begin
    shine_spectrum_specs = Symbol[]
    show_shine_Tb_spectrum && push!(shine_spectrum_specs, :brightness)
    show_shine_tau_spectrum && push!(shine_spectrum_specs, :optical_depth)
    if isempty(shine_spectrum_specs)
        fig_shine_spectrum = Figure(size = (900, 180))
        Label(fig_shine_spectrum[1, 1], L"\mathrm{Select\ at\ least\ one\ H\,I\ spectrum.}", fontsize = 20)
    else
        fig_shine_spectrum = Figure(size = (570length(shine_spectrum_specs), 390))
        for (index, spectrum_type) in enumerate(shine_spectrum_specs)
            if spectrum_type == :brightness
                ax = latex_axis(fig_shine_spectrum[1, index],
                    xlabel = L"v_{\mathrm{LOS}}\;[\mathrm{km\,s}^{-1}]", ylabel = L"T_B\;[\mathrm{K}]")
                show_shine_Tb_spectrum && lines!(ax, shine_velocity_axis,
                    @view(shine_Tb[Int(shine_sky_i), Int(shine_sky_j), :]);
                    color = MHD_COLORS[1], linewidth = 2.5, label = "T_B")
                axislegend(ax; position = :rt, framevisible = false)
            else
                ax = latex_axis(fig_shine_spectrum[1, index],
                    xlabel = L"v_{\mathrm{LOS}}\;[\mathrm{km\,s}^{-1}]", ylabel = L"\tau_{21}")
                lines!(ax, shine_velocity_axis,
                    @view(shine_tau[Int(shine_sky_i), Int(shine_sky_j), :]);
                    color = MHD_COLORS[5], linewidth = 2.5)
            end
        end
    end
    display_shine_spectrum ? fig_shine_spectrum : nothing
end

# ╔═╡ a0050005-6f8c-4d0c-9a10-000000000005
begin
    shine_structure_specs = [
        (data = shine_NHI, label = L"N_{\mathrm{HI}}\;[\mathrm{cm}^{-2}]", color = MHD_COLORS[1], period = nothing),
        (data = shine_NCNM, label = L"N_{\mathrm{CNM}}\;[\mathrm{cm}^{-2}]", color = MHD_COLORS[2], period = nothing),
        (data = shine_NLNM, label = L"N_{\mathrm{LNM}}\;[\mathrm{cm}^{-2}]", color = MHD_COLORS[3], period = nothing),
        (data = shine_NWNM, label = L"N_{\mathrm{WNM}}\;[\mathrm{cm}^{-2}]", color = MHD_COLORS[4], period = nothing),
        (data = shine_peakTb, label = L"\max_v T_B\;[\mathrm{K}]", color = MHD_COLORS[5], period = nothing),
        (data = shine_mom0, label = L"M_0\;[\mathrm{K\,km\,s}^{-1}]", color = MHD_COLORS[6], period = nothing),
        (data = shine_mom1, label = L"M_1\;[\mathrm{km\,s}^{-1}]", color = MHD_COLORS[1], period = nothing),
        (data = shine_mom2, label = L"M_2\;[\mathrm{km\,s}^{-1}]", color = MHD_COLORS[2], period = nothing),
        (data = shine_peak_tau, label = L"\max_v\tau_{21}", color = MHD_COLORS[3], period = nothing),
        (data = shine_fftcnm, label = L"f_{\mathrm{CNM}}^{\mathrm{FFT}}", color = MHD_COLORS[4], period = nothing),
        (data = shine_rgb_blue_map, label = L"W_{\mathrm{blue}}\;[\mathrm{K\,km\,s}^{-1}]", color = :dodgerblue3, period = nothing),
        (data = shine_rgb_green_map, label = L"W_{\mathrm{green}}\;[\mathrm{K\,km\,s}^{-1}]", color = :seagreen3, period = nothing),
        (data = shine_rgb_red_map, label = L"W_{\mathrm{red}}\;[\mathrm{K\,km\,s}^{-1}]", color = :firebrick3, period = nothing),
    ]
    fig_shine_structure = display_observational_structure_functions ?
        observational_structure_figure(shine_structure_specs, cube, sky_dims,
            observational_structure_order, observational_structure_samples;
            heading = "SHINE observable structure functions") : Figure(size = (900, 120))
    display_observational_structure_functions ? fig_shine_structure : nothing
end

# ╔═╡ e1000001-6f8c-4d0c-9a10-000000000001
begin
    export_figure_options = [
        "shine" => "SHINE H I maps",
        "shine_structure" => "SHINE observable structure functions",
        "shine_rgb" => "SHINE velocity RGB composite",
        "shine_spectrum" => "SHINE H I spectrum",
    ]
    export_figure_registry = Dict(
        "shine" => fig_shine,
        "shine_structure" => fig_shine_structure,
        "shine_rgb" => fig_shine_rgb,
        "shine_spectrum" => fig_shine_spectrum,
    )
    nothing
end

# ╔═╡ e1000002-6f8c-4d0c-9a10-000000000002
md"""
---

## Figure export

| Export setting | Control |
|:--|:--|
| Figure | $(@bind export_figure_key PlutoUI.Select(export_figure_options; default = "shine")) |
| Format | $(@bind export_figure_format PlutoUI.Select(["PNG", "PDF"]; default = "PDF")) |
"""

# ╔═╡ e1000003-6f8c-4d0c-9a10-000000000003
begin
    export_extension = lowercase(export_figure_format)
    export_mime = export_figure_format == "PNG" ? MIME"image/png"() : MIME"application/pdf"()
    export_buffer = IOBuffer()
    show(export_buffer, export_mime, export_figure_registry[export_figure_key])
    export_bytes = take!(export_buffer)
    export_run_slug = replace(lowercase(selected_run), r"[^a-z0-9]+" => "_")
    export_filename = "shine_$(export_figure_key)_$(export_run_slug)_snapshot_$(lpad(selected_snapshot, 3, '0')).$(export_extension)"
    PlutoUI.DownloadButton(export_bytes, export_filename)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FITSIO = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CairoMakie = "0.15"
FFTW = "1.10"
FITSIO = "0.17"
HDF5 = "0.17"
LaTeXStrings = "1.4"
Pluto = "1"
PlutoUI = "0.7"
StatsBase = "0.34"
julia = "1.11"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "d35da2752d645b5b15e5516bc40c59ff87b7b4e8"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "daa72978cd7a624246e894a4f4f067706d4e17e2"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.BitFlags]]
git-tree-sha1 = "bbe1079eecf9c9fbb52765193ad2bae27ae09bc8"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "8c6b984c3928736d455eb53a6adf881457825269"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.7.2"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d846411d1d7eb34739105f3d6fbba7c93beded3d"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.6.4+0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "47142129b1777e21da58cff265050b10d8560588"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.13"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1fa950ebc3e37eccd51c6a8fe1f92f7d86263522"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.7+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "eeaad7cef88554c2fa56b5a3f71cfd5cb708c662"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "7bc84b769c1d384315e7b5c4ac03a6c303e6cf35"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.8"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "4358750bb58a3caefd5f37a4a0c5bfdbbf075252"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.6"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "cd3c5ac74cd3923c8945c6a81518c46abd0e73a3"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.129"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e6c4a6407a949e79a9d3f249bf49e6987c80e01f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.2+0"

[[deps.ExpressionExplorer]]
git-tree-sha1 = "5f1c005ed214356bbe41d442cc1ccd416e510b7e"
uuid = "21656369-7473-754a-2065-74616d696c43"
version = "1.1.4"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "7a58e45171b63ed4782f2d36fdee8713a469e6e0"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.2+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "f57de3f533590c785210893030736dc11c4a4afb"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.17.5"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bad39456d9f0166184fce2248783dd9862645c1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Gamma]]
git-tree-sha1 = "86f86b6168a016ed88e4ae4e64577b98c3b59e8e"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.1.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "364685f5ffde25deb1bbcfd5bb278a5c6b7a9b37"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.11"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.GracefulPkg]]
deps = ["Compat", "Pkg", "TOML"]
git-tree-sha1 = "a854d6c0e9fb561b88cd20b4ad64f518cb1bfb8d"
uuid = "828d9ff0-206c-6161-646e-6576656f7244"
version = "2.4.3"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "491ea627ac824619f34168e29a0427a9e00e3e40"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.3"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "MPIABI_jll", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "aws_c_s3_jll", "dlfcn_win32_jll", "libaec_jll", "mpif_jll"]
git-tree-sha1 = "45337643a2d97262d5fe72ce1f13e8a662d13d62"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "2.1.2+0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "c35847ca5b4997fc8418836354a56c459bcf48d8"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.14.0+0"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "18d7deab5fb0440dc6a7b6993c5c27b25420de10"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.29"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "c3ee408ae340565f41699e3a3fa1053698c7626e"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.10"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1dae3057da6f2b9c857afef03177bbdc7c4afe92"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b7970cef8ae1c990ba0c09cd8bdc1145e006632f"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "22.1.7+0"

[[deps.LRUCache]]
git-tree-sha1 = "5519b95a490ff5fe629c4a7aa3b3dfc9160498b3"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.6.2"
weakdeps = ["Serialization"]

    [deps.LRUCache.extensions]
    SerializationExt = ["Serialization"]

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "0f2da712350b020bc3957f269c9caad516383ee0"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MPIABI_jll]]
deps = ["Artifacts", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9be143b6045719e8fb019d2b3bc2aebad1184fef"
uuid = "b5ada748-db0f-5fc0-8972-9331c762740c"
version = "0.1.5+0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "07dbec8aab01696edc0151a401a6cdfe95b9b885"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "5.0.1+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8e98d5d80b87403c311fd51e8455d4546ba7a5f8"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.12"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "675df097f8eeb28998b2cfe3b25655af73d5f7df"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.6+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "f2c8715d05bf10f9d4dc354e69dee30b6be53239"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.13"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.Malt]]
deps = ["Distributed", "Logging", "RelocatableFolders", "Serialization", "Sockets"]
git-tree-sha1 = "c2335b4e291f2422e2be8abf8936ccad58a98992"
uuid = "36869731-bdee-424d-aa32-cab38c994e3b"
version = "1.4.1"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e8dcbeef032ba2f9051a44ac22b4e54e3a1a0099"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.6"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dafdaa3ff15f20ff703d909d3a6f574a5b0586f3"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.33+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "0d621a4beb5e48d195f907c3c5b0bea285d9ff9d"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.13+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "6d6c0ca4824268c1a7dca1f4721c535ac63d9074"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.11+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d8cce34295c55f47be683580f44791716045b8fe"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.7+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "26766d4b5f1a410c218a19b85a672c6edb693c65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.40"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "32a4e09c5f29402573d673901778a0e03b0807b9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.6"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Downloads", "ExpressionExplorer", "FileWatching", "GracefulPkg", "HTTP", "HypertextLiteral", "InteractiveUtils", "LRUCache", "Logging", "LoggingExtras", "MIMEs", "Malt", "Markdown", "MsgPack", "Pkg", "PlutoDependencyExplorer", "PrecompileSignatures", "PrecompileTools", "REPL", "Random", "RegistryInstances", "RelocatableFolders", "SHA", "Scratch", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "fe7515cf6ddb62e738d924e4ca2dddaa60ff80ba"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "1.0.3"

[[deps.PlutoDependencyExplorer]]
deps = ["ExpressionExplorer", "InteractiveUtils", "Markdown"]
git-tree-sha1 = "c3e5073a977b1c58b2d55c1ec187c3737e64e6af"
uuid = "72656b73-756c-7461-726b-72656b6b696b"
version = "1.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "7fb25a964849d90a0446366cdefca822e0e84900"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.6"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "6547cbdd8ce32efba0d21c5a40fa96d1a3548f9f"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "770240df9a3b8888065046948f7a09b4e0f997d5"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "0f38a06c83f0007bbab3cf911262841c9a0f07e0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.13.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.aws_c_auth_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_cal_jll", "aws_c_http_jll", "aws_c_sdkutils_jll"]
git-tree-sha1 = "8cab83c96af80a1be968251ce1a0548a7545484d"
uuid = "2b3700d1-4306-52e2-a478-c162f0c514be"
version = "0.9.6+0"

[[deps.aws_c_cal_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "22c0f42f4a1f0dc5dcfa8fd267c4ac407c455e7a"
uuid = "70f11efc-bab2-57f1-b0f3-22aad4e67c4b"
version = "0.9.13+0"

[[deps.aws_c_common_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a759cb9bf456ad792cc7898a81ae333cce9ef02a"
uuid = "73048d1d-b8c4-5092-a58d-866c5e8d1e50"
version = "0.12.6+0"

[[deps.aws_c_compression_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "7910c72f45f44afd297c39fe43b99c56d5ed22ec"
uuid = "73a04cd5-f3d7-5bac-9290-e8adb709f224"
version = "0.3.2+0"

[[deps.aws_c_http_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_compression_jll", "aws_c_io_jll"]
git-tree-sha1 = "e358d5a001ef7afbd4f8c5225322512819cda2f2"
uuid = "3254fc65-9028-534d-aa9d-d76d128babc6"
version = "0.10.13+0"

[[deps.aws_c_io_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_cal_jll", "aws_c_common_jll", "s2n_tls_jll"]
git-tree-sha1 = "7e481d474b2087ee8bbf55b81bf9119f21e396d9"
uuid = "13c41daa-f319-5298-b5eb-5754e0170d52"
version = "0.26.3+0"

[[deps.aws_c_s3_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_auth_jll", "aws_c_common_jll", "aws_c_http_jll", "aws_checksums_jll", "s2n_tls_jll"]
git-tree-sha1 = "3e9917ab25114feba657e71be41cad068b9f6595"
uuid = "bd1f34fb-993f-5903-a121-aaf302eed6d4"
version = "0.11.5+0"

[[deps.aws_c_sdkutils_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "c43dfba2c1ab9ea9f02f2c80e86fa16f6460244e"
uuid = "1282aa60-004d-510b-9f52-12498d409daa"
version = "0.2.4+1"

[[deps.aws_checksums_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "aws_c_common_jll"]
git-tree-sha1 = "2570c8e23f4771a087b12a47edcaaa670ac05a01"
uuid = "b2a88e68-78e7-5e94-8c20-c02986ec140e"
version = "0.2.10+0"

[[deps.dlfcn_win32_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e141d67ffe550eadfb5af1bdbdaf138031e4805f"
uuid = "c4b69c83-5512-53e3-94e6-de98773c479f"
version = "1.4.2+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60f4792734488db6f42e2c7699f1d4594780bd03"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.7+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.mpif_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIABI_jll", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "TOML"]
git-tree-sha1 = "a8083ee0737c243c8f40a4ba86a0956997facb73"
uuid = "9aeb927a-4695-514f-a259-621a69f20ec0"
version = "0.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.s2n_tls_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "64ae051c6f03044eb7d98027d1b552b4e21e650c"
uuid = "cddc5d3d-934d-5d3a-9747-62fc12ea3f48"
version = "1.7.3+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─34d6b5f4-9c2d-42d5-9034-543aeb8ae151
# ╟─8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686
# ╠═3ef88702-bef5-4eca-a151-df97aa7ec2c4
# ╟─bdc11245-d76c-45fe-b79d-7b64861f5f53
# ╟─7bd6f2c9-ae49-4636-a251-f526ab347125
# ╟─55cbdbf4-e0f2-431e-a736-09f41ab7ee75
# ╟─98360288-85ca-4551-bdde-c12c7a329302
# ╟─6a9cfec2-2c80-4d72-94c0-cdb47aa4f046
# ╟─18e0d2d1-56ca-46cc-a13e-77f142612b5a
# ╟─e734297f-506e-45e1-8cb7-b2ae671893eb
# ╠═7174c31f-f186-48f9-b66e-29cf9a1c1fe3
# ╟─353cc6fb-c801-448a-a2c5-23dfd1541704
# ╟─a8ef96ab-0ddd-4eb2-a216-b7d96c2a9a08
# ╟─32110739-b60e-4592-856a-dd74f7a37401
# ╟─94a0a0dc-baf6-4e62-a51e-dc6124d98fd4
# ╟─c12d1f54-40b8-4865-9562-8dcb519f924a
# ╟─62440e86-b560-44ad-bb0a-43ae62e73fc3
# ╠═47b786d6-c7b5-44f4-946a-b8c485ad6380
# ╟─478ec2f3-e057-4720-809c-17ca0a3dac21
# ╟─bcc05889-02bb-47cf-b672-139e8efe4137
# ╠═5ad762e1-105f-4cf7-9cf1-e0bb8c6f1bf5
# ╟─27e51ba5-4592-4766-9dde-0de383a889a0
# ╠═a0050005-6f8c-4d0c-9a10-000000000005
# ╠═e1000001-6f8c-4d0c-9a10-000000000001
# ╠═e1000002-6f8c-4d0c-9a10-000000000002
# ╠═e1000003-6f8c-4d0c-9a10-000000000003
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
