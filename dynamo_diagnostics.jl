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

    "Return a memoized scalar summary, or `nothing` if it has not been computed."
    reduction_hit(key) = lock(() -> get(REDUCTION_CACHE, key, nothing), CACHE_LOCK)

    "Memoize a scalar summary and return it."
    store_reduction!(key, value) =
        lock(() -> (REDUCTION_CACHE[key] = value), CACHE_LOCK)

    """
    Return the cached raw cube for `key`, evaluating `build()` on a miss.

    The previous entry is evicted before a cache miss is read. Consequently a
    slider change never leaves several raw snapshots resident in this cache.
    """
    function cached_raw_cube!(key, build)
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
        lock(CACHE_LOCK) do
            empty!(RAW_CUBE_CACHE)
            empty!(RAW_CUBE_CACHE_BYTES)
            empty!(RAW_CUBE_CACHE_ORDER)
        end
        value = build()
        bytes = Base.summarysize(value)
        lock(CACHE_LOCK) do
            empty!(RAW_CUBE_CACHE)
            empty!(RAW_CUBE_CACHE_BYTES)
            empty!(RAW_CUBE_CACHE_ORDER)
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
# Interactive MHD diagnostics

Select the data, simulation, snapshot, and line of sight below. Run a result cell with **Shift+Enter**; Pluto evaluates its required dependencies automatically.
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
        half .* (lower .+ upper)
    end

    """
    Cheap fingerprint that changes whenever a snapshot is rewritten on disk.

    Including it in the cache keys means the caches invalidate themselves when a
    simulation is re-run, instead of serving stale arrays.
    """
    function snapshot_fingerprint(path)
        isdir(path) || return (mtime(path), filesize(path))
        Tuple((basename(entry), mtime(entry), filesize(entry))
            for entry in sort(readdir(path; join = true)))
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
        stems = Set(normalize_field_name(splitext(file)[1]) for file in readdir(path)
            if is_fits_file(joinpath(path, file)))
        any(alias -> normalize_field_name(alias) in stems, FIELD_ALIASES[:rho]) &&
            all(field -> any(alias -> normalize_field_name(alias) in stems,
                FIELD_ALIASES[field]), (:vx, :vy, :vz))
    end

    function snapshot_sources(cube_directory)
        isdir(cube_directory) || return String[]
        fits_directory_is_snapshot(cube_directory) && return [cube_directory]
        sources = String[]
        for path in readdir(cube_directory; join = true)
            if is_hdf5_file(path) || is_fits_file(path) ||
                    fits_directory_is_snapshot(path)
                push!(sources, path)
            end
        end
        sort(sources)
    end

    function expand_home(path)
        path == "~" && return homedir()
        startswith(path, "~/") ? joinpath(homedir(), path[3:end]) : path
    end

    function discover_cube_directories(path; depth = 0, max_depth = 32)
        isdir(path) || return String[]
        direct_snapshots = snapshot_sources(path)
        isempty(direct_snapshots) || return [abspath(path)]
        depth >= max_depth && return String[]
        found = String[]
        for entry in sort(readdir(path))
            startswith(entry, ".") && continue
            child = joinpath(path, entry)
            isdir(child) || continue
            append!(found, discover_cube_directories(child;
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

# ╔═╡ 3df9ad9a-b865-41f5-8d7a-34a52d0292dd
Markdown.parse("""
### Repository status

| Item | Active value |
|:--|:--|
| Comparison parameter | $(comparison_parameter) |
| Resolved data root | **$(ROOT)** |
| Snapshot limit | **$(MAX_SNAPSHOTS_PER_RUN)** per simulation, evenly spaced |
| Discovered runs | $(run_summary) |
""")

# ╔═╡ 290ecb0a-880a-44bd-9ddc-e6ee12b41a06
md"""
### Navigation and physical units

**Run** selects the active cube. **Simulations in comparative plots** starts with that single run and can contain any number of simulations.
"""

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
    active_time_value = snapshot_time(
        run_files[selected_run][selected_snapshot]) * time_unit_Myr
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
    function read_raw_cube(path)
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
        raw_fields, raw_length, raw_time = with_hdf5_file(path) do h
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
    load_raw_cube(path) = cached_raw_cube!(raw_cube_key(path), () -> read_raw_cube(path))

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

# ╔═╡ 5fbd0d53-09e8-4b39-bb1d-29c8cc45c6ee
begin
    validation_fields = [
        "Density ρ" => cube.rho,
        "Pressure P" => cube.P,
        "Velocity vx" => cube.vx,
        "Velocity vy" => cube.vy,
        "Velocity vz" => cube.vz,
        "Magnetic field Bx" => cube.bx,
        "Magnetic field By" => cube.by,
        "Magnetic field Bz" => cube.bz,
    ]
    cube_cell_count = length(cube.rho)
    validation_rows = String[]
    total_invalid_values = 0
    for (field_name, values) in validation_fields
        finite_count = count(isfinite, values)
        invalid_count = length(values) - finite_count
        total_invalid_values += invalid_count
        finite_fraction = 100finite_count / max(length(values), 1)
        finite_fraction_text = @sprintf("%.3f", finite_fraction)
        push!(validation_rows,
            "| $(field_name) | $(finite_fraction_text)% | $(invalid_count) |")
    end
    invalid_density_count = count(value -> !isfinite(value) || value <= 0, cube.rho)
    invalid_pressure_count = count(value -> !isfinite(value) || value <= 0, cube.P)
    cube_memory_mib = sum(Base.summarysize(last(pair)) for pair in validation_fields) / 2.0^20
    validation_state = total_invalid_values == 0 && invalid_density_count == 0 &&
        invalid_pressure_count == 0 ? "All required fields are finite and physically positive." :
        "Invalid cells are masked by finite-only statistics; inspect the counts below before interpretation."
    source_text = replace(String(selected_path), "`" => "\\`")
    Markdown.parse("""
    ### Active cube validation

    **$(validation_state)**

    | Property | Value |
    |:--|:--|
    | Source | **$(source_text)** |
    | Grid shape | **$(join(size(cube.rho), " × "))** |
    | Physical box | **$(@sprintf("%.5g × %.5g × %.5g", cube.L...))** ``\\mathrm{pc}^3`` |
    | Approximate field memory | **$(@sprintf("%.2f", cube_memory_mib)) MiB** |
    | Non-positive or non-finite density cells | **$(invalid_density_count)** / **$(cube_cell_count)** |
    | Non-positive or non-finite pressure cells | **$(invalid_pressure_count)** / **$(cube_cell_count)** |

    | Required field | Finite values | Invalid values |
    |:--|--:|--:|
    $(join(validation_rows, "\n"))
    """)
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

# ╔═╡ 76d249f9-d6fd-4513-aa76-7fb386058c37
begin
    heatmap_specs = NamedTuple[]
    projected_B_lic = show_projected_B_lic ? lic_texture(Bsky1, Bsky2;
        niter = Int(lic_iterations), len = Int(lic_length),
        normalize_vectors = lic_normalize_vectors,
        amplitude_weight = lic_amplitude_weight,
        amplitude_floor = Float64(lic_amplitude_floor), seed = Int(lic_seed)) : nothing

    function add_heatmap!(enabled, A, label, colormap, use_log; signed = false, overlay_B = false)
        enabled || return
        transformed = use_log ? (signed ? symlog10(A) : safe_log10.(A)) : A
        latex_label = use_log ?
            (signed ? latexstring("\\operatorname{symlog}_{10}\\left(", label, "\\right)") :
                latexstring("\\log_{10}\\left(", label, "\\right)")) :
            latexstring(label)
        push!(heatmap_specs, (
            data = transformed,
            label = latex_label,
            colormap = signed ? :balance : colormap,
            diverging = signed,
            overlay_B = overlay_B,
        ))
    end

    add_heatmap!(show_column, column_density, "N_{\\mathrm{H}}/(\\mathrm{cm}^{-2})", :magma, log_column;
        overlay_B = show_projected_B || show_projected_B_lic)
    add_heatmap!(show_temperature, Tmean, "T/\\mathrm{K}", :thermal, log_temperature)
    add_heatmap!(show_blos, Blos, "B_{\\mathrm{LOS}}/\\mu\\mathrm{G}", :balance, log_blos; signed = true)
    add_heatmap!(show_bsky, Bsky, "B_{\\mathrm{sky}}/\\mu\\mathrm{G}", :viridis, log_bsky)
    add_heatmap!(show_Bmean, Bmean_projected, "\\langle |B|\\rangle/\\mu\\mathrm{G}", :viridis, log_Bmean)
    add_heatmap!(show_velocity, velocity_projected, "\\langle |\\delta v|\\rangle/(\\mathrm{km\\,s}^{-1})", :plasma, log_velocity)
    add_heatmap!(show_vorticity, omega_map, "\\langle |\\omega|\\rangle/\\mathrm{Myr}^{-1}", :inferno, log_vorticity)
    add_heatmap!(show_magnetic_energy, magnetic_energy_map,
        "\\Sigma_{E,\\mathrm{mag}}/(\\mathrm{erg\\,cm}^{-2})", :magma, log_magnetic_energy)
    add_heatmap!(show_kinetic_energy, kinetic_energy_map,
        "\\Sigma_{E,\\mathrm{kin}}/(\\mathrm{erg\\,cm}^{-2})", :viridis, log_kinetic_energy)
    add_heatmap!(show_thermal_energy, thermal_energy_map,
        "\\Sigma_{E,\\mathrm{therm}}/(\\mathrm{erg\\,cm}^{-2})", :thermal, log_thermal_energy)

    if isempty(heatmap_specs)
        fig_maps = Figure(size = (900, 180))
        Label(fig_maps[1, 1], L"\mathrm{Select\ at\ least\ one\ field\ to\ display\ a\ heatmap.}", fontsize = 20)
    else
        ncols = length(heatmap_specs) == 1 ? 1 : 2
        nrows = cld(length(heatmap_specs), ncols)
        fig_maps = Figure(size = (560ncols, 420nrows + 60))
        for (index, spec) in enumerate(heatmap_specs)
            row, col = cld(index, ncols), mod1(index, ncols)
            panel = fig_maps[row, col] = GridLayout()
            ax = latex_axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            colorrange = robust_colorrange(spec.data, color_percentile; diverging = spec.diverging)
            hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap, colorrange)
            if spec.overlay_B && show_projected_B_lic
                lic_colors = [RGBAf(0, 0, 0, 0),
                    RGBAf(1, 1, 1, Float64(lic_opacity))]
                heatmap!(ax, sky_coordinates[1], sky_coordinates[2], projected_B_lic;
                    colormap = lic_colors, colorrange = (0.0, 1.0), transparency = true)
            end
            latex_colorbar(panel[1, 2], hm, label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 20)

            if spec.overlay_B && show_projected_B
                stride = arrow_stride
                ix = collect(1:stride:size(Bsky1, 1))
                iy = collect(1:stride:size(Bsky1, 2))
                points = [Point2f(sky_coordinates[1][i], sky_coordinates[2][j]) for i in ix for j in iy]
                dirs = [Vec2f(Bsky1[i, j], Bsky2[i, j]) for i in ix for j in iy]
                arrow_length_pc = 0.65stride * step(sky_coordinates[1])
                arrows2d!(ax, points, dirs; normalize = true, lengthscale = arrow_length_pc,
                    align = :center, color = (:white, 0.85), shaftwidth = 1.5,
                    tipwidth = 7, tiplength = 5)
            end
        end
    end
    display_projected_maps ? fig_maps : nothing
end

# ╔═╡ 72626861-42ce-4ac0-b980-78f498f8a629
md"""
---

## 3. Physical probability density functions

The panels compare every run selected in **Simulations in comparative plots** at the selected snapshot index (clamped to the last available snapshot of shorter runs). They show number density $n$, magnetic-field strength $|B|$, and turbulent speed $|\delta\mathbf v|$ in physical units. Shared $\log_{10}X$ bins are used across runs. Their vertical axes are probability densities per dex and satisfy $\int (\mathrm{d}\mathcal P/\mathrm{d}\log_{10}X)\,\mathrm{d}\log_{10}X=1$.

**Display physical PDFs:** $(@bind display_pdfs PlutoUI.CheckBox(default = true))

| PDF panel | Display |
|:--|:--:|
| Number-density PDF | $(@bind show_pdf_density PlutoUI.CheckBox(default = true)) |
| Magnetic-field PDF | $(@bind show_pdf_magnetic PlutoUI.CheckBox(default = true)) |
| Turbulent-speed PDF | $(@bind show_pdf_velocity PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 496cbf2d-77a1-4a1a-b760-d4f8ea2ea9de
begin
    function density_pdf(values, weights, edges)
        valid = isfinite.(values) .& isfinite.(weights) .& (weights .>= 0)
        values, weights = Float64.(values[valid]), Float64.(weights[valid])
        isempty(values) && return (Float64[], Float64[])
        h = fit(Histogram, values, Weights(weights), edges)
        centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
        pdf = h.weights ./ max(sum(h.weights .* diff(edges)), eps())
        centers, pdf
    end

    function cube_pdf_samples(c)
        local_mag = magnetic_fields(c)
        local_turb = turbulent_velocity(c)
        local_weights = pdf_weighting == "mass" ?
            vec(Float64.(c.rho)) : ones(length(c.rho))
        local_B = GAUSS_TO_MICROGAUSS .* local_mag.B
        local_v = sqrt.(local_turb.dv2)
        local_mean_B = finite_mean(local_mag.B)
        (
            density = vec(safe_log10.(number_density(c.rho))),
            magnetic = vec(safe_log10.(local_B)),
            velocity = vec(safe_log10.(local_v)),
            normalized_B = vec(safe_log10.(local_mag.B ./ local_mean_B)),
            weights = local_weights,
        )
    end

    function comparative_pdf(samples_by_run, field, bins)
        combined = vcat([finite_values(getfield(samples_by_run[label], field))
            for label in comparison_run_labels]...)
        isempty(combined) && return Dict{String, Any}()
        lo, hi = quantile(combined, (0.001, 0.999))
        lo == hi && ((lo, hi) = (lo - 0.5, hi + 0.5))
        edges = range(lo, hi; length = Int(bins) + 1)
        Dict(label => density_pdf(getfield(samples_by_run[label], field),
                samples_by_run[label].weights, edges)
            for label in comparison_run_labels)
    end

    number_density_cells = number_density(cube.rho)
    magnetic_strength_uG = GAUSS_TO_MICROGAUSS .* mag.B
    turbulent_speed_kms = sqrt.(turb.dv2)
    mean_B = finite_mean(mag.B)
    sB = vec(safe_log10.(mag.B ./ mean_B))
    logn = vec(safe_log10.(number_density_cells))
    logBphysical = vec(safe_log10.(magnetic_strength_uG))
    logvphysical = vec(safe_log10.(turbulent_speed_kms))
    comparative_pdf_samples = Dict(label => cube_pdf_samples(comparison_cube(label))
        for label in comparison_run_labels)
    density_pdfs = comparative_pdf(comparative_pdf_samples, :density, nbins)
    magnetic_pdfs = comparative_pdf(comparative_pdf_samples, :magnetic, nbins)
    velocity_pdfs = comparative_pdf(comparative_pdf_samples, :velocity, nbins)
    normalized_B_pdfs = comparative_pdf(comparative_pdf_samples, :normalized_B, nbins)
end

# ╔═╡ 1e0e6c0e-1ae0-40a1-a6f8-9c18fae91961
begin
    pdf_specs = NamedTuple[]
    show_pdf_density && push!(pdf_specs, (pdfs = density_pdfs,
        xlabel = L"n\;[\mathrm{cm}^{-3}]"))
    show_pdf_magnetic && push!(pdf_specs, (pdfs = magnetic_pdfs,
        xlabel = L"|B|\;[\mu\mathrm{G}]"))
    show_pdf_velocity && push!(pdf_specs, (pdfs = velocity_pdfs,
        xlabel = L"|\delta v|\;[\mathrm{km\,s}^{-1}]"))
    if isempty(pdf_specs)
        fig_pdf = Figure(size = (900, 180))
        Label(fig_pdf[1, 1], L"\mathrm{Select\ at\ least\ one\ PDF.}", fontsize = 20)
    else
        fig_pdf = Figure(size = (420length(pdf_specs), 430))
        for (j, spec) in enumerate(pdf_specs)
            ax = latex_axis(fig_pdf[1, j], xlabel = spec.xlabel,
                ylabel = L"\mathrm{d}\mathcal{P}/\mathrm{d}\log_{10}X", xscale = log10,
                xticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), xminorticksvisible = true)
            for label in comparison_run_labels
                logx, probability = spec.pdfs[label]
                stairs!(ax, 10.0 .^ logx, probability;
                    color = run_colors[label], linewidth = 2.5, step = :center,
                    label = legend_run_label(label))
            end
            ylims!(ax, low = 0)
        end
        Legend(fig_pdf[2, 1:length(pdf_specs)],
            [LineElement(color = run_colors[label], linewidth = 2.5)
                for label in comparison_run_labels],
            legend_run_label.(comparison_run_labels), L"\mathrm{Simulation}";
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    display_pdfs ? fig_pdf : nothing
end

# ╔═╡ fa155a62-da75-4530-b5e9-215fd4f66412
md"""
---

## 4. Thermodynamic phase diagram

This comparative figure shows one joint distribution of number density $n$ and thermal pressure $P/k_{\mathrm B}$ for each selected simulation. The plotted coordinates are $\log_{10}n$ and $\log_{10}(P/k_{\mathrm B})$; empty probability bins are masked and all panels share one color scale.

The conversion follows the units defined in section 1: $n=\rho/(\mu m_{\mathrm H})$ and $P/k_{\mathrm B}$ in $\mathrm{K\,cm}^{-3}$. The selected **PDF weighting** is applied consistently. The optional Koyama--Inutsuka equilibrium curve satisfies $n\Lambda(T)=\Gamma$ with $P/k_{\mathrm B}=nT$.

**Display the pressure-density phase diagram:** $(@bind display_phase_diagram PlutoUI.CheckBox(default = true))

| Phase-diagram setting | Control |
|:--|:--|
| Bins per axis | $(@bind phase_bins PlutoUI.Slider(20:5:120; default = 60, show_value = true)) |
| Thermal-equilibrium curve | $(@bind show_phase_equilibrium PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 9da572fe-21f2-43df-9320-b8742fd64773
begin
    function phase_histogram(logn, logpk, weights, bins)
        valid = isfinite.(logn) .& isfinite.(logpk) .& isfinite.(weights) .& (weights .>= 0)
        x, y, w = logn[valid], logpk[valid], weights[valid]
        isempty(x) && error("No finite positive cells are available for the phase diagram.")
        xlo, xhi = quantile(x, (0.001, 0.999))
        ylo, yhi = quantile(y, (0.001, 0.999))
        xlo == xhi && ((xlo, xhi) = (xlo - 0.5, xhi + 0.5))
        ylo == yhi && ((ylo, yhi) = (ylo - 0.5, yhi + 0.5))
        xedges = range(xlo, xhi; length = bins + 1)
        yedges = range(ylo, yhi; length = bins + 1)
        histogram = zeros(Float64, bins, bins)
        for index in eachindex(x)
            ix = searchsortedlast(xedges, x[index])
            iy = searchsortedlast(yedges, y[index])
            1 <= ix <= bins && 1 <= iy <= bins || continue
            histogram[ix, iy] += w[index]
        end
        dx, dy = step(xedges), step(yedges)
        probability = histogram ./ max(sum(histogram) * dx * dy, eps())
        log_probability = map(value -> value > 0 ? log10(value) : NaN, probability)
        xcenters = (xedges[1:end-1] .+ xedges[2:end]) ./ 2
        ycenters = (yedges[1:end-1] .+ yedges[2:end]) ./ 2
        (; xcenters, ycenters, log_probability)
    end

    function koyama_inutsuka_equilibrium(; temperature_min = 10.0,
            temperature_max = 1.0e5, samples = 2400)
        temperature_K = 10.0 .^ range(log10(temperature_min), log10(temperature_max);
            length = samples)
        cooling_over_heating_cm3 =
            1.0e7 .* exp.(-1.184e5 ./ (temperature_K .+ 1000.0)) .+
            1.4e-2 .* sqrt.(temperature_K) .* exp.(-92.0 ./ temperature_K)
        equilibrium_density_cm3 = 1.0 ./ cooling_over_heating_cm3
        equilibrium_pressure_over_k = equilibrium_density_cm3 .* temperature_K
        valid = isfinite.(equilibrium_density_cm3) .&
            isfinite.(equilibrium_pressure_over_k) .&
            (equilibrium_density_cm3 .> 0) .& (equilibrium_pressure_over_k .> 0)
        (
            logn = safe_log10.(equilibrium_density_cm3[valid]),
            logpk = safe_log10.(equilibrium_pressure_over_k[valid]),
        )
    end

    phase_data_by_run = Dict(label => begin
        local_cube = comparison_cube(label)
        local_logn = vec(safe_log10.(number_density(local_cube.rho)))
        local_logpk = vec(safe_log10.(local_cube.P ./ K_B_CGS))
        local_weights = pdf_weighting == "mass" ?
            vec(Float64.(local_cube.rho)) : ones(length(local_cube.rho))
        phase_histogram(local_logn, local_logpk, local_weights, phase_bins)
    end for label in comparison_run_labels)
    phase_equilibrium = koyama_inutsuka_equilibrium()
end

# ╔═╡ 41b4eb12-889d-43b3-87c2-fc7cccf8679f
begin
    phase_panel_count = length(comparison_run_labels)
    fig_phase = Figure(size = (560phase_panel_count + 80, 520))
    combined_phase_probability = vcat([
        finite_values(phase_data_by_run[label].log_probability)
        for label in comparison_run_labels]...)
    phase_range = robust_colorrange(combined_phase_probability, 99.0)
    phase_heatmap = nothing
    for (panel_index, label) in enumerate(comparison_run_labels)
        phase_data = phase_data_by_run[label]
        phase_axis = latex_axis(fig_phase[1, panel_index],
            xlabel = L"\log_{10}\!\left(n/\mathrm{cm}^{-3}\right)",
            ylabel = L"\log_{10}\!\left[(P/k_B)/(\mathrm{K\,cm}^{-3})\right]",
            title = latex_run_label(label))
        phase_heatmap = heatmap!(phase_axis, phase_data.xcenters, phase_data.ycenters,
            phase_data.log_probability; colormap = :magma, colorrange = phase_range)
        if show_phase_equilibrium
            lines!(phase_axis, phase_equilibrium.logn, phase_equilibrium.logpk;
                color = :black, linewidth = 4.5)
            lines!(phase_axis, phase_equilibrium.logn, phase_equilibrium.logpk;
                color = :white, linewidth = 2.5,
                label = L"n\Lambda(T)=\Gamma")
            panel_index == 1 && axislegend(phase_axis; position = :rt, labelsize = 14)
        end
        phase_dx = length(phase_data.xcenters) > 1 ?
            phase_data.xcenters[2] - phase_data.xcenters[1] : 1.0
        phase_dy = length(phase_data.ycenters) > 1 ?
            phase_data.ycenters[2] - phase_data.ycenters[1] : 1.0
        xlims!(phase_axis, first(phase_data.xcenters) - phase_dx / 2,
            last(phase_data.xcenters) + phase_dx / 2)
        ylims!(phase_axis, first(phase_data.ycenters) - phase_dy / 2,
            last(phase_data.ycenters) + phase_dy / 2)
    end
    latex_colorbar(fig_phase[1, phase_panel_count + 1], phase_heatmap,
        label = L"\log_{10}\mathcal{P}_{2\mathrm{D}}", tickformat = latex_ticklabels)
    colsize!(fig_phase.layout, phase_panel_count + 1, 22)
    display_phase_diagram ? fig_phase : nothing
end

# ╔═╡ 71c8ea26-d2ad-4430-9265-0b28d45bba1c
md"""
### Interval magnetic-growth rate

The local magnetic-growth rate is evaluated between consecutive snapshots rather than with one fit over the complete selected window:

```math
\Gamma_{B,i+1/2}
=\frac{\ln B_{i+1}-\ln B_i}{t_{i+1}-t_i}.
```

Time, sonic Mach number, and Alfvénic Mach number are assigned their midpoint values over the same interval. Positive $\Gamma_B$ indicates magnetic amplification; negative $\Gamma_B$ indicates decay.

**Display interval growth-rate relations:** $(@bind display_gamma_relations PlutoUI.CheckBox(default = false))

| Growth-rate setting | Control |
|:--|:--|
| Magnetic field | $(@bind gamma_relation_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Simulations | Global multi-selection in section 1 |
| $\Gamma_B(t)$ | $(@bind show_gamma_time PlutoUI.CheckBox(default = true)) |
| $\Gamma_B(\mathcal{M})$ | $(@bind show_gamma_mach PlutoUI.CheckBox(default = true)) |
| $\Gamma_B(\mathcal{M}_{\mathrm A})$ | $(@bind show_gamma_alfven_mach PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ cd7e037c-62f5-4682-92c2-92af7169d692
md"""
### Multi-snapshot normalized magnetic evolution

Selected snapshots are superimposed in the same figure using

```math
\ln\!\left(\frac{B_i}{B_0}\right),
```

where $B_0$ is the first available snapshot of each run. Color identifies the run; every snapshot uses the same circular marker to keep temporal comparisons readable. All axes remain linear because the logarithm is applied to the magnetic-field ratio itself.

**Display multi-snapshot normalized magnetic evolution:** $(@bind display_normalized_B_relations PlutoUI.CheckBox(default = false))

| Normalized-field relation | Control |
|:--|:--|
| Magnetic field | $(@bind normalized_B_relation_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Snapshots displayed | $(@bind normalized_B_snapshot_indices PlutoUI.MultiSelect(collect(1:maximum_snapshot_count); default = unique([1, cld(maximum_snapshot_count, 2), maximum_snapshot_count]))) |
| Simulations displayed | Global multi-selection in section 1 |
| $\ln(B/B_0)$ versus $t$ | $(@bind show_normalized_B_time PlutoUI.CheckBox(default = true)) |
| $\ln(B/B_0)$ versus $\mathcal{M}$ | $(@bind show_normalized_B_mach PlutoUI.CheckBox(default = true)) |
| $\ln(B/B_0)$ versus $\mathcal{M}_{\mathrm A}$ | $(@bind show_normalized_B_alfven_mach PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ dcc8f8f3-daaf-4d4b-92a3-4919ed5e36de
md"""
---

## 6. Global time evolution

Selected panels compare every discovered run as a function of physical time. Turbulent velocity is measured after subtracting the mass-weighted bulk motion.

$\mathcal{M}=v_{\mathrm{rms}}/c_{s,\mathrm{rms}}$ measures compressibility, while $\mathcal{M}_{\mathrm A}=v_{\mathrm{rms}}/v_{\mathrm A,\mathrm{rms}}$ compares turbulence with Alfvén-wave propagation. Values above unity are supersonic or super-Alfvénic, respectively. Dotted magnetic-field curves show the theoretical exponentials selected in section 5.

**Display global time evolution:** $(@bind display_global_evolution PlutoUI.CheckBox(default = false))

| Time-evolution panel | Display |
|:--|:--:|
| Sonic Mach number | $(@bind show_time_mach PlutoUI.CheckBox(default = true)) |
| Alfvénic Mach number | $(@bind show_time_alfven PlutoUI.CheckBox(default = true)) |
| Magnetic field | $(@bind show_time_magnetic PlutoUI.CheckBox(default = true)) |
| Magnetic-to-kinetic energy ratio | $(@bind show_time_energy_ratio PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 62bb58f1-37c3-4adb-a7fc-939f71a56635
md"""
### Magnetic field by thermal phase

This diagnostic follows the magnetic-field strength after dividing every snapshot into three temperature phases. The phase boundaries are adjustable. **Volume average** assigns the same weight to every cell, while **Density-weighted average** uses $w_i=\rho_i$. The displayed phase mean and RMS are

```math
\langle |B|\rangle_{w,\phi}=\frac{\sum_{i\in\phi}w_i|B_i|}{\sum_{i\in\phi}w_i},\qquad
B_{\mathrm{rms},w,\phi}=\left(\frac{\sum_{i\in\phi}w_iB_i^2}{\sum_{i\in\phi}w_i}\right)^{1/2},
```

with $w_i=1$ for volume normalization or $w_i=\rho_i$ for density normalization. Values are shown in $\mu\mathrm G$; no additional normalization by a global $B_0$ is applied in this phase diagnostic.

```math
\mathrm{CNM}:\ T<T_{\mathrm C},\qquad
\mathrm{LNM}:\ T_{\mathrm C}\leq T<T_{\mathrm W},\qquad
\mathrm{WNM}:\ T\geq T_{\mathrm W}.
```

**Display magnetic field by phase:** $(@bind display_phase_B_time PlutoUI.CheckBox(default = false))

| Phase-field setting | Control |
|:--|:--|
| Cold/Lukewarm boundary $T_{\mathrm C}$ [$\mathrm{K}$] | $(@bind phase_cold_boundary_K PlutoUI.NumberField(10.0:10.0:5000.0; default = 200.0)) |
| Lukewarm/Warm boundary $T_{\mathrm W}$ [$\mathrm{K}$] | $(@bind phase_warm_boundary_K PlutoUI.NumberField(100.0:50.0:50000.0; default = 2000.0)) |
| Magnetic statistic | $(@bind phase_B_statistic PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Phase normalization | $(@bind phase_B_weighting PlutoUI.Select(["Volume average" => "volume", "Density-weighted average" => "density"]; default = "volume")) |
| Logarithmic $B$ axis | $(@bind log_phase_B_time PlutoUI.CheckBox(default = true)) |
| CNM curve | $(@bind show_phase_B_cold PlutoUI.CheckBox(default = true)) |
| LNM curve | $(@bind show_phase_B_lukewarm PlutoUI.CheckBox(default = true)) |
| WNM curve | $(@bind show_phase_B_warm PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 904ba663-d536-4b27-a379-4af927b0affb
begin
    """
    Bulk dynamics of one cube, accumulated in a single pass over the grid.

    The broadcast form of this reduction allocated roughly thirty cube-sized
    temporaries; the loop allocates none.
    """
    function bulk_metrics_from_cube(c, gamma)
        m = magnetic_fields(c)
        u = turbulent_velocity(c)
        rho, P, B, B2, dv2 = c.rho, c.P, m.B, m.B2, u.dv2
        inverse_kms2 = 1 / KM_CM^2
        thermal_factor = gamma > 1 + sqrt(eps(Float64)) ? 1 / (gamma - 1) : 1.0
        mass, mass_dv2 = 0.0, 0.0
        cs_mass, cs_weighted = 0.0, 0.0
        va_mass, va_weighted = 0.0, 0.0
        Ekin, Emag, Etherm = 0.0, 0.0, 0.0
        B_sum, B_count = 0.0, 0
        B2_sum, B2_count = 0.0, 0
        @inbounds for index in eachindex(rho)
            density = Float64(rho[index])
            pressure = Float64(P[index])
            magnitude = Float64(B[index])
            square = Float64(B2[index])
            turbulence = Float64(dv2[index])

            if isfinite(magnitude)
                B_sum += magnitude
                B_count += 1
            end
            if isfinite(square)
                B2_sum += square
                B2_count += 1
            end

            magnetic_energy = square / (8pi)
            isfinite(magnetic_energy) && (Emag += magnetic_energy)
            kinetic_energy = 0.5 * density * turbulence * KM_CM^2
            isfinite(kinetic_energy) && (Ekin += kinetic_energy)
            thermal_energy = thermal_factor * pressure
            isfinite(thermal_energy) && (Etherm += thermal_energy)

            (isfinite(density) && density > 0 && isfinite(turbulence)) || continue
            mass += density
            mass_dv2 += density * turbulence
            cs2 = gamma * pressure / density * inverse_kms2
            if isfinite(cs2) && cs2 >= 0
                cs_mass += density
                cs_weighted += density * cs2
            end
            va2 = square / (4pi * density) * inverse_kms2
            if isfinite(va2) && va2 >= 0
                va_mass += density
                va_weighted += density * va2
            end
        end
        vrms = mass > 0 ? sqrt(mass_dv2 / mass) : NaN
        cs_rms = cs_mass > 0 ? sqrt(cs_weighted / cs_mass) : NaN
        va_rms = va_mass > 0 ? sqrt(va_weighted / va_mass) : NaN
        mean_B = B_count > 0 ? B_sum / B_count : NaN
        mean_B2 = B2_count > 0 ? B2_sum / B2_count : NaN
        (
            t = c.t,
            mach = vrms / max(cs_rms, eps()),
            mach_alfven = vrms / max(va_rms, eps()),
            Bmean = GAUSS_TO_MICROGAUSS * mean_B,
            Brms = GAUSS_TO_MICROGAUSS * sqrt(mean_B2),
            energy_ratio = Ekin > 0 ? Emag / Ekin : NaN,
            kin_mag = Ekin > 0 && Emag > 0 ? Ekin / Emag : NaN,
            therm_mag = Etherm > 0 && Emag > 0 ? Etherm / Emag : NaN,
            kin_therm = Ekin > 0 && Etherm > 0 ? Ekin / Etherm : NaN,
        )
    end

    """
    Mean and RMS field of the three thermal phases, in a single pass over the cube.

    Matches the mask-based form this replaces, including its treatment of
    non-finite temperatures: `NaN` fails all three comparisons and so belongs to
    no phase, while `-Inf` falls in the cold phase and `+Inf` in the warm one.
    """
    function phase_metrics_from_cube(c, molecular_weight, cold_boundary_K,
            warm_boundary_K, phase_weighting)
        m = magnetic_fields(c)
        rho, P, B, B2 = c.rho, c.P, m.B, m.B2
        cold_boundary = min(cold_boundary_K, warm_boundary_K)
        warm_boundary = max(cold_boundary_K, warm_boundary_K)
        temperature_factor = molecular_weight * M_H_CGS / K_B_CGS
        density_weighted = phase_weighting == "density"
        counts = zeros(Int, 3)
        weights = zeros(3)
        weighted_B = zeros(3)
        weighted_B2 = zeros(3)
        @inbounds for index in eachindex(rho)
            density = Float64(rho[index])
            magnitude = Float64(B[index])
            square = Float64(B2[index])
            (isfinite(density) && isfinite(magnitude) && isfinite(square)) || continue
            temperature = temperature_factor * Float64(P[index]) / density
            isnan(temperature) && continue
            phase = temperature < cold_boundary ? 1 :
                temperature < warm_boundary ? 2 : 3
            counts[phase] += 1
            density_weighted && density <= 0 && continue
            weight = density_weighted ? density : 1.0
            weights[phase] += weight
            weighted_B[phase] += weight * magnitude
            weighted_B2[phase] += weight * square
        end
        statistics = ntuple(3) do phase
            (counts[phase] == 0 || weights[phase] <= 0) && return (mean = NaN, rms = NaN)
            (
                mean = GAUSS_TO_MICROGAUSS * weighted_B[phase] / weights[phase],
                rms = GAUSS_TO_MICROGAUSS * sqrt(weighted_B2[phase] / weights[phase]),
            )
        end
        (
            t = c.t,
            B_cold_mean = statistics[1].mean,
            B_cold_rms = statistics[1].rms,
            B_lukewarm_mean = statistics[2].mean,
            B_lukewarm_rms = statistics[2].rms,
            B_warm_mean = statistics[3].mean,
            B_warm_rms = statistics[3].rms,
        )
    end

    """
    Sweep one run once, deriving both series from a single read of each snapshot.

    Each summary is memoized on its own parameters, so changing `gamma` leaves
    the phase series untouched and changing a phase boundary leaves the bulk
    series untouched. Only the missing half is recomputed, and the cube behind it
    is read only when neither half is already cached.
    """
    function run_metric_series(paths, gamma, molecular_weight, cold_boundary_K,
            warm_boundary_K, phase_weighting, with_phase)
        bulk = Vector{Any}(undef, length(paths))
        phase = with_phase ? Vector{Any}(undef, length(paths)) : nothing
        keys_by_index = map(paths) do path
            signature = cube_signature(path)
            ((:bulk, signature, gamma),
             (:phase, signature, molecular_weight, cold_boundary_K,
                warm_boundary_K, phase_weighting))
        end

        # Serve everything already memoized first, so that only the snapshots
        # that still have to be read reach the parallel section below.
        pending = Int[]
        for index in eachindex(paths)
            bulk_key, phase_key = keys_by_index[index]
            bulk_hit = reduction_hit(bulk_key)
            phase_hit = with_phase ? reduction_hit(phase_key) : nothing
            if isnothing(bulk_hit) || (with_phase && isnothing(phase_hit))
                push!(pending, index)
            else
                bulk[index] = bulk_hit
                with_phase && (phase[index] = phase_hit)
            end
        end

        # read_raw_cube bypasses the single-entry cube cache on purpose: a sweep
        # visits each snapshot once, so caching here would only evict the
        # interactively selected cube, and would have parallel workers evict
        # each other.
        function reduce_snapshot!(index)
            bulk_key, phase_key = keys_by_index[index]
            raw = read_raw_cube(paths[index])
            c = scale_raw_cube(raw)
            bulk[index] = store_reduction!(bulk_key, bulk_metrics_from_cube(c, gamma))
            with_phase && (phase[index] = store_reduction!(phase_key,
                phase_metrics_from_cube(c, molecular_weight, cold_boundary_K,
                    warm_boundary_K, phase_weighting)))
            Base.summarysize(raw)
        end

        if !isempty(pending)
            # The first snapshot is read serially, both to size the remaining
            # workers against free memory and to surface a malformed file as a
            # plain error rather than one wrapped in a task failure.
            cube_bytes = reduce_snapshot!(first(pending))
            rest = @view pending[2:end]
            isempty(rest) ||
                parallel_foreach(reduce_snapshot!, rest, sweep_concurrency(cube_bytes))
        end
        (bulk = identity.(bulk), phase = with_phase ? identity.(phase) : nothing)
    end

    temporal_series_requested = display_global_evolution ||
        display_gamma_relations || display_normalized_B_relations ||
        display_growth_fit || display_phase_B_time || display_energy_time ||
        display_polarization_time

    function metadata_only_series(label)
        [(
            t = Float64(time_unit_Myr) * Float64(time),
            mach = NaN, mach_alfven = NaN, Bmean = NaN, Brms = NaN,
            energy_ratio = NaN, kin_mag = NaN, therm_mag = NaN,
            kin_therm = NaN,
        ) for time in run_times[label]]
    end

    metric_series_by_run = temporal_series_requested ? Dict(
            label => run_metric_series(run_files[label], Float64(gamma),
                Float64(mean_molecular_weight), Float64(phase_cold_boundary_K),
                Float64(phase_warm_boundary_K), String(phase_B_weighting),
                display_phase_B_time && label in comparison_run_labels)
            for label in analysis_series_labels
        ) : Dict(
            label => (bulk = metadata_only_series(label), phase = nothing)
            for label in analysis_series_labels
        )

    all_series = Dict(label => metric_series_by_run[label].bulk
        for label in analysis_series_labels)

    phase_B_series_by_run = display_phase_B_time ? Dict(
        label => metric_series_by_run[label].phase
        for label in comparison_run_labels
    ) : Dict{String, Any}()
end

# ╔═╡ 36aef377-3de7-435a-af83-3a90421e3159
md"""
---

## 5. Magnetic amplification and growth-rate fit

The fitted physical model is $B(t)=A\exp[\Gamma_B(t-t_0)]$, where $t_0$ is the first valid snapshot time and both $A$ and $\Gamma_B$ are inferred from the snapshots in the **Fit window**. Equivalently, the regression is $\ln B=a+\Gamma_B(t-t_0)$ with a free intercept. The shaded interval marks the selected fit range. This avoids treating the first magnetic-field measurement as exact or giving it special leverage in the regression.

The figures show $\ln(B/B_0)$ on linear axes, with $B_0$ the first valid measured field used only as a plotting normalization. The reported $R^2$ uses the usual mean-centred total sum of squares appropriate to this intercept fit. The theoretical comparison curves remain anchored at the first valid measured snapshot. Because $E_{\mathrm B}\propto B^2$, magnetic energy grows at rate $2\Gamma_B$.

**Display the magnetic growth-rate fit:** $(@bind display_growth_fit PlutoUI.CheckBox(default = false))

| Setting | Control |
|:--|:--|
| Fitted field | $(@bind growth_fit_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Fit window (snapshot indices) | $(@bind growth_fit_window PlutoUI.RangeSlider(1:length(run_files[selected_run]); default = min(2, length(run_files[selected_run])):min(4, length(run_files[selected_run])), show_value = true)) |
| Theoretical $\Gamma_{B,1}$ [$\mathrm{Myr}^{-1}$] | $(@bind theory_gamma_1 PlutoUI.NumberField(-0.50:0.001:0.50; default = 0.01)) |
| Theoretical $\Gamma_{B,2}$ [$\mathrm{Myr}^{-1}$] | $(@bind theory_gamma_2 PlutoUI.NumberField(-0.50:0.001:0.50; default = 0.03)) |
| Theoretical $\Gamma_{B,3}$ [$\mathrm{Myr}^{-1}$] | $(@bind theory_gamma_3 PlutoUI.NumberField(-0.50:0.001:0.50; default = 0.05)) |
| Data and fitted exponential panel | $(@bind show_growth_fit_panel PlutoUI.CheckBox(default = true)) |
| Theoretical-growth comparison panel | $(@bind show_growth_theory_panel PlutoUI.CheckBox(default = true)) |
| Display fitted curve | $(@bind show_growth_fit PlutoUI.CheckBox(default = true)) |
| Logarithmic $B(t)$ axis in section 6 | $(@bind log_B_time PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 89c33295-34e7-49ec-8d04-52b2aac29cff
begin
    growth_series = all_series[selected_run]
    growth_times = Float64.(getfield.(growth_series, :t))
    growth_B = growth_fit_field == "Mean field ⟨B⟩" ?
        Float64.(getfield.(growth_series, :Bmean)) :
        Float64.(getfield.(growth_series, :Brms))
    growth_reference_index = findfirst(i -> isfinite(growth_times[i]) &&
        isfinite(growth_B[i]) && growth_B[i] > 0, eachindex(growth_B))
    growth_B0 = isnothing(growth_reference_index) ? NaN : growth_B[growth_reference_index]
    growth_t0 = isnothing(growth_reference_index) ? NaN : growth_times[growth_reference_index]
    growth_log_B = map(growth_B) do value
        isfinite(value) && value > 0 ? log(value) : NaN
    end
    growth_elapsed_time = growth_times .- growth_t0
    growth_indices = collect(growth_fit_window)
    growth_indices = filter(i -> 1 <= i <= length(growth_times) &&
        isfinite(growth_times[i]) && isfinite(growth_elapsed_time[i]) &&
        isfinite(growth_log_B[i]), growth_indices)

    function exponential_growth_fit(elapsed_time, log_field, indices)
        if length(indices) < 2
            return (log_amplitude = NaN, gamma = NaN, gamma_error = NaN, r2 = NaN)
        end
        xfit, yfit = elapsed_time[indices], log_field[indices]
        xmean, ymean = mean(xfit), mean(yfit)
        centered_x = xfit .- xmean
        sxx = sum(abs2, centered_x)
        if sxx <= 0
            return (log_amplitude = NaN, gamma = NaN, gamma_error = NaN, r2 = NaN)
        end
        gamma_fit = sum(centered_x .* (yfit .- ymean)) / sxx
        log_amplitude = ymean - gamma_fit * xmean
        residuals = yfit .- (log_amplitude .+ gamma_fit .* xfit)
        ssres = sum(abs2, residuals)
        sstot = sum(abs2, yfit .- mean(yfit))
        r2 = sstot > 0 ? 1 - ssres / sstot : NaN
        gamma_error = length(indices) > 2 ?
            sqrt((ssres / (length(indices) - 2)) / sxx) : NaN
        (; log_amplitude, gamma = gamma_fit, gamma_error, r2)
    end

    growth_fit = exponential_growth_fit(growth_elapsed_time, growth_log_B, growth_indices)
    theory_gammas = Float64[theory_gamma_1, theory_gamma_2, theory_gamma_3]
    theory_B_curves = [growth_B0 .* exp.(Γ .* (growth_times .- growth_t0)) for Γ in theory_gammas]
    fitted_B_curve = exp.(growth_fit.log_amplitude .+ growth_fit.gamma .* growth_elapsed_time)
    fitted_ratio_curve = fitted_B_curve ./ growth_B0
    growth_has_interval = !isempty(growth_indices)
    growth_first_index = growth_has_interval ? first(growth_indices) : missing
    growth_last_index = growth_has_interval ? last(growth_indices) : missing
    growth_gamma_text = string(round(growth_fit.gamma; sigdigits = 5))
    growth_error_text = string(round(growth_fit.gamma_error; sigdigits = 3))
    growth_r2_text = string(round(growth_fit.r2; sigdigits = 4))
    growth_amplitude_text = string(round(exp(growth_fit.log_amplitude); sigdigits = 5))
    growth_energy_gamma_text = string(round(2 * growth_fit.gamma; sigdigits = 5))
end

# ╔═╡ 8120f7e9-74ed-4a48-b2f4-dbc75ebf0132
Markdown.parse("""
### Current fit result

| Quantity | Value |
|:--|:--|
| Run and field | **$(selected_run)** — $(growth_fit_field) |
| Fitted snapshot interval | **$(growth_first_index)–$(growth_last_index)** |
| Fitted amplitude ``A`` at ``t_0`` | **$(growth_amplitude_text)** ``\\mu\\mathrm{G}`` |
| Magnetic growth rate ``\\Gamma_B`` | **$(growth_gamma_text)** ``\\mathrm{Myr}^{-1}`` |
| Standard error ``\\sigma_{\\Gamma_B}`` | **$(growth_error_text)** ``\\mathrm{Myr}^{-1}`` |
| Coefficient of determination ``R^2`` | **$(growth_r2_text)** |
| Magnetic-energy growth rate ``2\\Gamma_B`` | **$(growth_energy_gamma_text)** ``\\mathrm{Myr}^{-1}`` |
""")

# ╔═╡ 3bed2efb-f849-4ddc-9f49-ed5e3d683370
begin
    theory_colors = MHD_COLORS[4:6]
    time_panel_keys = Symbol[]
    show_time_mach && push!(time_panel_keys, :mach)
    show_time_alfven && push!(time_panel_keys, :alfven)
    show_time_magnetic && push!(time_panel_keys, :magnetic)
    show_time_energy_ratio && push!(time_panel_keys, :energy)
    if isempty(time_panel_keys)
        fig_time = Figure(size = (900, 180))
        Label(fig_time[1, 1], L"\mathrm{Select\ at\ least\ one\ time-evolution\ panel.}", fontsize = 20)
    else
        time_ncols = length(time_panel_keys) == 1 ? 1 : 2
        time_nrows = cld(length(time_panel_keys), time_ncols)
        fig_time = Figure(size = (560time_ncols, 380time_nrows + 110))
        time_axes = Dict{Symbol, Any}()
        for (index, key) in enumerate(time_panel_keys)
            row, col = cld(index, time_ncols), mod1(index, time_ncols)
            time_axes[key] = key == :mach ?
                latex_axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]", ylabel = L"\mathcal{M}") :
                key == :alfven ?
                latex_axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]", ylabel = L"\mathcal{M}_{\mathrm{A}}") :
                key == :magnetic ?
                latex_axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]",
                    ylabel = L"B\;[\mu\mathrm{G}]", yscale = log_B_time ? log10 : identity,
                    yticks = log_B_time ? DECADE_TICKS : Makie.automatic,
                    yminorticks = log_B_time ? IntervalsBetween(9) : IntervalsBetween(5),
                    yminorticksvisible = true) :
                latex_axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]",
                    ylabel = L"E_B/E_{\mathrm{kin}}", yscale = log10,
                    yticks = DECADE_TICKS,
                    yminorticks = IntervalsBetween(9), yminorticksvisible = true)
        end
        for label in comparison_run_labels
            s = all_series[label]
            ts = getfield.(s, :t)
            if haskey(time_axes, :mach)
                lines!(time_axes[:mach], ts, getfield.(s, :mach);
                    color = run_colors[label], linewidth = 2)
                scatter!(time_axes[:mach], ts, getfield.(s, :mach);
                    color = run_colors[label], markersize = 7)
            end
            haskey(time_axes, :alfven) &&
                lines!(time_axes[:alfven], ts, getfield.(s, :mach_alfven);
                    color = run_colors[label], linewidth = 2)
            if haskey(time_axes, :magnetic)
                lines!(time_axes[:magnetic], ts, getfield.(s, :Bmean);
                    color = run_colors[label], linewidth = 2)
                lines!(time_axes[:magnetic], ts, getfield.(s, :Brms);
                    color = run_colors[label], linestyle = :dash)
            end
            haskey(time_axes, :energy) &&
                lines!(time_axes[:energy], ts, getfield.(s, :energy_ratio);
                    color = run_colors[label], linewidth = 2)
        end
        if haskey(time_axes, :magnetic)
            for (Γ, curve, color) in zip(theory_gammas, theory_B_curves, theory_colors)
                lines!(time_axes[:magnetic], growth_times, curve;
                    color, linewidth = 2, linestyle = :dot)
            end
            if show_growth_fit && isfinite(growth_fit.gamma)
                lines!(time_axes[:magnetic], growth_times, fitted_B_curve;
                    color = :black, linewidth = 2.5, linestyle = :dashdot)
            end
            growth_has_interval && vspan!(time_axes[:magnetic],
                growth_times[growth_first_index], growth_times[growth_last_index];
                color = (:gray50, 0.10))
            ylims!(time_axes[:magnetic], high = 15.0)
        end

        time_legend_layout = GridLayout(fig_time[time_nrows + 1, 1:time_ncols])
        time_legend_column = 1
        Legend(time_legend_layout[1, time_legend_column],
            [LineElement(color = run_colors[label], linewidth = 2.5)
                for label in comparison_run_labels],
            legend_run_label.(comparison_run_labels), "Simulation";
            orientation = :horizontal, tellheight = true, framevisible = false,
            labelsize = 11)
        if haskey(time_axes, :magnetic)
            time_legend_column += 1
            Legend(time_legend_layout[1, time_legend_column], [
                LineElement(color = :gray25, linewidth = 2.5, linestyle = :solid),
                LineElement(color = :gray25, linewidth = 2.5, linestyle = :dash),
            ], ["⟨B⟩", "Bᵣₘₛ"], "Magnetic statistic";
                orientation = :horizontal, tellheight = true, framevisible = false,
                labelsize = 11)

            model_elements = Any[
                LineElement(color = color, linewidth = 2.2, linestyle = :dot)
                for color in theory_colors
            ]
            model_labels = Any[legend_rate_label(Γ) for Γ in theory_gammas]
            if show_growth_fit && isfinite(growth_fit.gamma)
                push!(model_elements,
                    LineElement(color = :black, linewidth = 2.5, linestyle = :dashdot))
                push!(model_labels, legend_rate_label(growth_fit.gamma; fitted = true))
            end
            time_legend_column += 1
            Legend(time_legend_layout[1, time_legend_column],
                model_elements, model_labels, "Selected-run model";
                orientation = :horizontal, tellheight = true, framevisible = false,
                labelsize = 11)
        end
    end
    stable_pluto_figure(display_global_evolution, fig_time)
end

# ╔═╡ f29d5dff-8810-4935-b880-9951bc1239fd
begin
    function interval_growth_rate(series, field)
        times = Float64.(getfield.(series, :t))
        magnetic_field = Float64.(getfield.(series, field))
        mach = Float64.(getfield.(series, :mach))
        mach_alfven = Float64.(getfield.(series, :mach_alfven))
        length(times) < 2 && return (
            time = Float64[], gamma = Float64[], mach = Float64[], mach_alfven = Float64[])
        delta_time = diff(times)
        gamma_interval = diff(log.(magnetic_field)) ./ delta_time
        interval_time = (times[1:end-1] .+ times[2:end]) ./ 2
        interval_mach = (mach[1:end-1] .+ mach[2:end]) ./ 2
        interval_mach_alfven = (mach_alfven[1:end-1] .+ mach_alfven[2:end]) ./ 2
        valid = isfinite.(gamma_interval) .& isfinite.(interval_time) .&
            isfinite.(interval_mach) .& isfinite.(interval_mach_alfven) .&
            (delta_time .> 0) .& (magnetic_field[1:end-1] .> 0) .&
            (magnetic_field[2:end] .> 0)
        (
            time = interval_time[valid],
            gamma = gamma_interval[valid],
            mach = interval_mach[valid],
            mach_alfven = interval_mach_alfven[valid],
        )
    end

    gamma_panel_specs = NamedTuple[]
    show_gamma_time && push!(gamma_panel_specs,
        (field = :time, xlabel = L"t\;[\mathrm{Myr}]"))
    show_gamma_mach && push!(gamma_panel_specs,
        (field = :mach, xlabel = L"\mathcal{M}"))
    show_gamma_alfven_mach && push!(gamma_panel_specs,
        (field = :mach_alfven, xlabel = L"\mathcal{M}_{\mathrm A}"))
    gamma_runs = comparison_run_labels
    gamma_field_symbol = gamma_relation_field == "Mean field ⟨B⟩" ? :Bmean : :Brms

    if isempty(gamma_panel_specs)
        fig_gamma_relations = Figure(size = (900, 180))
        Label(fig_gamma_relations[1, 1],
            L"\mathrm{Select\ at\ least\ one\ growth-rate\ relation.}", fontsize = 20)
    else
        gamma_ncols = length(gamma_panel_specs) == 1 ? 1 : min(3, length(gamma_panel_specs))
        fig_gamma_relations = Figure(size = (430gamma_ncols, 450))
        for (panel_index, spec) in enumerate(gamma_panel_specs)
            gamma_axis = latex_axis(fig_gamma_relations[1, panel_index],
                xlabel = spec.xlabel, ylabel = L"\Gamma_B\;[\mathrm{Myr}^{-1}]")
            hlines!(gamma_axis, [0.0]; color = (:gray35, 0.65),
                linestyle = :dash, linewidth = 1.5)
            for label in gamma_runs
                relation = interval_growth_rate(all_series[label], gamma_field_symbol)
                horizontal = getfield(relation, spec.field)
                lines!(gamma_axis, horizontal, relation.gamma;
                    color = run_colors[label], linewidth = 2.2, label = legend_run_label(label))
                scatter!(gamma_axis, horizontal, relation.gamma;
                    color = run_colors[label], markersize = 7)
            end
        end
        Legend(fig_gamma_relations[2, 1:length(gamma_panel_specs)],
            [LineElement(color = run_colors[label], linewidth = 2.4) for label in gamma_runs],
            legend_run_label.(gamma_runs); orientation = :horizontal,
            tellheight = true, framevisible = false)
    end
    display_gamma_relations ? fig_gamma_relations : nothing
end

# ╔═╡ 7a4472c3-48c5-44c7-a487-d1021b84ee3a
begin
    normalized_B_panel_specs = NamedTuple[]
    show_normalized_B_time && push!(normalized_B_panel_specs,
        (field = :t, xlabel = L"t\;[\mathrm{Myr}]"))
    show_normalized_B_mach && push!(normalized_B_panel_specs,
        (field = :mach, xlabel = L"\mathcal{M}"))
    show_normalized_B_alfven_mach && push!(normalized_B_panel_specs,
        (field = :mach_alfven, xlabel = L"\mathcal{M}_{\mathrm A}"))
    normalized_B_runs = comparison_run_labels
    normalized_B_field_symbol = normalized_B_relation_field == "Mean field ⟨B⟩" ? :Bmean : :Brms
    normalized_B_valid_snapshots = sort(unique(filter(index -> index >= 1,
        Int.(normalized_B_snapshot_indices))))

    if isempty(normalized_B_panel_specs) || isempty(normalized_B_valid_snapshots) ||
            isempty(normalized_B_runs)
        fig_normalized_B_relations = Figure(size = (900, 180))
        Label(fig_normalized_B_relations[1, 1],
            L"\mathrm{Select\ at\ least\ one\ relation\ and\ one\ snapshot.}", fontsize = 20)
    else
        normalized_B_ncols = length(normalized_B_panel_specs)
        fig_normalized_B_relations = Figure(size = (430normalized_B_ncols, 500))
        for (panel_index, spec) in enumerate(normalized_B_panel_specs)
            axis = latex_axis(fig_normalized_B_relations[1, panel_index],
                xlabel = spec.xlabel, ylabel = L"\ln(B/B_0)")
            hlines!(axis, [0.0]; color = (:gray35, 0.65), linestyle = :dash, linewidth = 1.4)
            for label in normalized_B_runs
                series = all_series[label]
                magnetic_field = Float64.(getfield.(series, normalized_B_field_symbol))
                B0 = first(magnetic_field)
                available = filter(index -> index <= length(series) &&
                    isfinite(magnetic_field[index]) && magnetic_field[index] > 0 &&
                    isfinite(B0) && B0 > 0, normalized_B_valid_snapshots)
                isempty(available) && continue
                horizontal_all = Float64.(getfield.(series, spec.field))
                horizontal = horizontal_all[available]
                logarithmic_ratio = log.(magnetic_field[available] ./ B0)
                valid = isfinite.(horizontal) .& isfinite.(logarithmic_ratio)
                lines!(axis, horizontal[valid], logarithmic_ratio[valid];
                    color = run_colors[label], linewidth = 2.2)
                scatter!(axis, horizontal[valid], logarithmic_ratio[valid];
                    color = run_colors[label], marker = :circle, markersize = 7)
            end
        end
        run_elements = [LineElement(color = run_colors[label], linewidth = 2.4)
            for label in normalized_B_runs]
        Legend(fig_normalized_B_relations[2, 1:length(normalized_B_panel_specs)],
            run_elements, legend_run_label.(normalized_B_runs), L"\mathrm{Simulation}";
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    display_normalized_B_relations ? fig_normalized_B_relations : nothing
end

# ╔═╡ a7def784-6c52-4f98-9d76-ece082b45229
begin
    normalized_B = growth_B ./ growth_B0
    ln_normalized_B = log.(max.(normalized_B, floatmin(Float64)))
    growth_panel_count = count(identity, (show_growth_fit_panel, show_growth_theory_panel))
    if growth_panel_count == 0
        fig_growth = Figure(size = (900, 180))
        Label(fig_growth[1, 1], L"\mathrm{Select\ at\ least\ one\ magnetic-growth\ panel.}", fontsize = 20)
    else
        fig_growth = Figure(size = (550growth_panel_count, 520))
        growth_column = 0
        if show_growth_fit_panel
            growth_column += 1
            ag1 = latex_axis(fig_growth[1, growth_column], xlabel = L"t\;[\mathrm{Myr}]",
                ylabel = L"\ln(B/B_0)")
            lines!(ag1, growth_times, ln_normalized_B; color = run_colors[selected_run],
                linewidth = 2.5)
            scatter!(ag1, growth_times, ln_normalized_B; color = run_colors[selected_run], markersize = 8)
            fit_legend_elements = Any[
                LineElement(color = run_colors[selected_run], linewidth = 2.5),
            ]
            fit_legend_labels = LaTeXString[L"\mathrm{Data}"]
            if show_growth_fit && isfinite(growth_fit.gamma)
                fit_t = growth_times[growth_indices]
                fit_log_ratio = growth_fit.log_amplitude - log(growth_B0) .+
                    growth_fit.gamma .* growth_elapsed_time[growth_indices]
                lines!(ag1, fit_t, fit_log_ratio; color = :black, linewidth = 3,
                    linestyle = :dash)
                push!(fit_legend_elements,
                    LineElement(color = :black, linewidth = 3, linestyle = :dash))
                push!(fit_legend_labels, latexstring(
                    "\\mathrm{Fit}:\\;\\Gamma_B=",
                    @sprintf("%.4g", growth_fit.gamma),
                    "\\;\\mathrm{Myr}^{-1}"))
            end
            growth_has_interval && vspan!(ag1,
                growth_times[growth_first_index], growth_times[growth_last_index];
                color = (:gray50, 0.10))
            Legend(fig_growth[2, growth_column], fit_legend_elements,
                fit_legend_labels; orientation = :horizontal, nbanks = 1,
                tellheight = true, framevisible = false, labelsize = 14)
        end
        if show_growth_theory_panel
            growth_column += 1
            ag2 = latex_axis(fig_growth[1, growth_column], xlabel = L"t\;[\mathrm{Myr}]",
                ylabel = L"\ln(B/B_0)")
            lines!(ag2, growth_times, ln_normalized_B; color = run_colors[selected_run],
                linewidth = 2.5)
            scatter!(ag2, growth_times, ln_normalized_B; color = run_colors[selected_run], markersize = 7)
            theory_legend_elements = Any[
                LineElement(color = run_colors[selected_run], linewidth = 2.5),
            ]
            theory_legend_labels = LaTeXString[L"\mathrm{Data}"]
            for (Γ, curve, color) in zip(theory_gammas, theory_B_curves, theory_colors)
                lines!(ag2, growth_times, log.(curve ./ growth_B0); color, linewidth = 2,
                    linestyle = :dot)
                push!(theory_legend_elements,
                    LineElement(color = color, linewidth = 2, linestyle = :dot))
                push!(theory_legend_labels, latexstring(
                    "\\Gamma_B=", @sprintf("%.3g", Γ),
                    "\\;\\mathrm{Myr}^{-1}"))
            end
            if show_growth_fit && isfinite(growth_fit.gamma)
                lines!(ag2, growth_times, log.(fitted_ratio_curve); color = :black,
                    linewidth = 2.5, linestyle = :dashdot)
                push!(theory_legend_elements,
                    LineElement(color = :black, linewidth = 2.5, linestyle = :dashdot))
                push!(theory_legend_labels, latexstring(
                    "\\mathrm{Best\\ fit}:\\;\\Gamma_B=",
                    @sprintf("%.4g", growth_fit.gamma),
                    "\\;\\mathrm{Myr}^{-1}"))
            end
            growth_has_interval && vspan!(ag2,
                growth_times[growth_first_index], growth_times[growth_last_index];
                color = (:gray50, 0.10))
            Legend(fig_growth[2, growth_column], theory_legend_elements,
                theory_legend_labels; orientation = :horizontal, nbanks = 2,
                tellheight = true, framevisible = false, labelsize = 13)
        end
        rowgap!(fig_growth.layout, 8)
    end
    display_growth_fit ? fig_growth : nothing
end

# ╔═╡ d56a8ba3-aa42-4351-a065-87275032a342
begin
    phase_B_specs = NamedTuple[]
    show_phase_B_cold && push!(phase_B_specs,
        (key = :cold, label = "CNM", linestyle = :solid))
    show_phase_B_lukewarm && push!(phase_B_specs,
        (key = :lukewarm, label = "LNM", linestyle = :dash))
    show_phase_B_warm && push!(phase_B_specs,
        (key = :warm, label = "WNM", linestyle = :dot))

    if !display_phase_B_time
        fig_phase_B_time = Figure(size = (900, 180))
    elseif isempty(phase_B_specs)
        fig_phase_B_time = Figure(size = (900, 180))
        Label(fig_phase_B_time[1, 1],
            L"\mathrm{Select\ at\ least\ one\ thermal\ phase.}", fontsize = 20)
    else
        fig_phase_B_time = Figure(size = (760, 480))
        phase_B_ylabel = phase_B_statistic == "Mean field ⟨B⟩" ?
            L"\langle |B|\rangle_{\mathrm{phase}}\;[\mu\mathrm{G}]" :
            L"B_{\mathrm{rms,phase}}\;[\mu\mathrm{G}]"
        phase_B_axis = latex_axis(fig_phase_B_time[1, 1],
            xlabel = L"t\;[\mathrm{Myr}]", ylabel = phase_B_ylabel,
            yscale = log_phase_B_time ? log10 : identity,
            yticks = log_phase_B_time ? DECADE_TICKS : Makie.automatic,
            yminorticks = log_phase_B_time ? IntervalsBetween(9) : IntervalsBetween(5),
            yminorticksvisible = true)
        statistic_suffix = phase_B_statistic == "Mean field ⟨B⟩" ? "mean" : "rms"
        for label in comparison_run_labels, spec in phase_B_specs
            phase_B_series = phase_B_series_by_run[label]
            phase_B_times = Float64.(getfield.(phase_B_series, :t))
            field = Symbol("B_", spec.key, "_", statistic_suffix)
            values = Float64.(getfield.(phase_B_series, field))
            valid = isfinite.(phase_B_times) .& isfinite.(values) .&
                (log_phase_B_time ? values .> 0 : trues(length(values)))
            lines!(phase_B_axis, phase_B_times[valid], values[valid];
                color = run_colors[label], linestyle = spec.linestyle, linewidth = 2.6)
            scatter!(phase_B_axis, phase_B_times[valid], values[valid];
                color = run_colors[label], marker = :circle, markersize = 7)
        end
        Legend(fig_phase_B_time[2, 1],
            [[LineElement(color = run_colors[label], linewidth = 2.5)
                for label in comparison_run_labels],
             [LineElement(color = :gray25, linestyle = spec.linestyle, linewidth = 2.5)
                for spec in phase_B_specs]],
            [legend_run_label.(comparison_run_labels), [spec.label for spec in phase_B_specs]],
            ["Simulation", "Phase"];
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    display_phase_B_time ? fig_phase_B_time : nothing
end

# ╔═╡ 1be45ca2-bd2f-4b77-9188-5b338d41483b
md"""
---

## 7. Normalized magnetic-field distribution

The map shows the active run's line-of-sight mean of $\log_{10}(B/\langle B\rangle)$. The PDF compares every selected simulation using shared bins and every cell of the corresponding cube. A value of zero corresponds to each cube's mean field strength; positive and negative values identify locally stronger and weaker fields.

**Display the normalized magnetic-field distribution:** $(@bind display_normalized_field PlutoUI.CheckBox(default = true))

| Normalized-field panel | Display |
|:--|:--:|
| $\log_{10}(B/\langle B\rangle)$ map | $(@bind show_logB_map PlutoUI.CheckBox(default = true)) |
| $\log_{10}(B/\langle B\rangle)$ PDF | $(@bind show_logB_histogram PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 89420b8d-c72e-4a04-91dd-043bc9ecef2e
begin
    logB = safe_log10.(mag.B ./ finite_mean(mag.B))
    logB_map = finite_mean_dims(logB, los_dim)
    logB_panel_count = count(identity, (show_logB_map, show_logB_histogram))
    if logB_panel_count == 0
        fig_logB = Figure(size = (900, 180))
        Label(fig_logB[1, 1], L"\mathrm{Select\ a\ normalized-field\ panel.}", fontsize = 20)
    else
        fig_logB = Figure(size = (530logB_panel_count, 390))
        logB_column = 0
        if show_logB_map
            logB_column += 1
            map_panel = fig_logB[1, logB_column] = GridLayout()
            axmap = latex_axis(map_panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            limB = max(maximum(abs, filter(isfinite, vec(logB_map)); init = 0.0),
                sqrt(eps(Float64)))
            hb = heatmap!(axmap, sky_coordinates[1], sky_coordinates[2], logB_map;
                colormap = :balance, colorrange = (-limB, limB))
            latex_colorbar(map_panel[1, 2], hb,
                label = L"\log_{10}(B/\langle B\rangle)", tickformat = latex_ticklabels)
            colsize!(map_panel, 2, 22)
        end
        if show_logB_histogram
            logB_column += 1
            axhist = latex_axis(fig_logB[1, logB_column],
                xlabel = L"\log_{10}(B/\langle B\rangle)", ylabel = L"\mathcal{P}")
            for label in comparison_run_labels
                logB_ratio_x, logB_ratio_p = normalized_B_pdfs[label]
                lines!(axhist, logB_ratio_x, logB_ratio_p;
                    color = run_colors[label], linewidth = 2.5,
                    label = legend_run_label(label))
            end
            axislegend(axhist; position = :rt, framevisible = false)
        end
    end
    display_normalized_field ? fig_logB : nothing
end

# ╔═╡ b1000001-6f8c-4d0c-9a10-000000000001
md"""
---

## Magnetic field--density analysis

| Diagnostic | Control |
|:--|:--|
| Display the comparative $B$--$n$ relation | $(@bind display_bn_relation PlutoUI.CheckBox(default = true)) |
| Minimum density used by the $B\propto n^\kappa$ fit [$\mathrm{cm}^{-3}$; $0$ uses all bins] | $(@bind bn_fit_min_density PlutoUI.NumberField(default = 0.0)) |
| Display the 3D HRO | $(@bind display_hro PlutoUI.CheckBox(default = true)) |
| Number of HRO density intervals | $(@bind hro_density_bin_count PlutoUI.Slider(4:1:20; default = 10, show_value = true)) |
| Display the projected HOG | $(@bind display_hog PlutoUI.CheckBox(default = true)) |
| HOG Gaussian smoothing FWHM [pixels] | $(@bind hog_smoothing_fwhm_pix PlutoUI.NumberField(0.0:0.5:20.0; default = 2.0)) |
| HOG gradient rejection percentile | $(@bind hog_gradient_percentile PlutoUI.Slider(0.0:5.0:90.0; default = 20.0, show_value = true)) |
| Apply HOG to $\log_{10}$ maps | $(@bind hog_logarithmic_maps PlutoUI.CheckBox(default = true)) |

HRO compares $\mathbf B$ with three-dimensional iso-density structures. HOG compares the gradients of projected $N_{\rm H}$ and density-weighted $|B|$. All selected simulations use the active snapshot index and line of sight.
"""

# ╔═╡ b1000002-6f8c-4d0c-9a10-000000000002
begin
    function magnetic_density_samples(c)
        local_n = number_density(c.rho)
        local_B = GAUSS_TO_MICROGAUSS .* magnetic_fields(c).B
        valid = isfinite.(local_n) .& isfinite.(local_B) .&
            (local_n .> 0) .& (local_B .> 0)
        (logn = log10.(Float64.(local_n[valid])),
            logB = log10.(Float64.(local_B[valid])))
    end

    function binned_magnetic_density(samples, logn_edges)
        centers = (logn_edges[1:end-1] .+ logn_edges[2:end]) ./ 2
        medians = fill(NaN, length(centers))
        lower = fill(NaN, length(centers))
        upper = fill(NaN, length(centers))
        counts = zeros(Int, length(centers))
        for bin in eachindex(centers)
            members = (samples.logn .>= logn_edges[bin]) .&
                (samples.logn .< logn_edges[bin + 1])
            values = samples.logB[members]
            counts[bin] = length(values)
            length(values) >= 4 || continue
            lower[bin], medians[bin], upper[bin] = quantile(values, (0.16, 0.50, 0.84))
        end
        (; centers, medians, lower, upper, counts)
    end

    function magnetic_density_fit(profile, minimum_density)
        valid = isfinite.(profile.centers) .& isfinite.(profile.medians) .&
            (profile.counts .>= 4)
        minimum_density > 0 &&
            (valid .&= profile.centers .>= log10(Float64(minimum_density)))
        x, y = profile.centers[valid], profile.medians[valid]
        length(x) >= 2 || return (slope = NaN, intercept = NaN)
        xmean, ymean = mean(x), mean(y)
        denominator = sum(abs2, x .- xmean)
        denominator > 0 || return (slope = NaN, intercept = NaN)
        slope = sum((x .- xmean) .* (y .- ymean)) / denominator
        (slope = slope, intercept = ymean - slope * xmean)
    end

    bn_samples_by_run = Dict(label => magnetic_density_samples(comparison_cube(label))
        for label in comparison_run_labels)
    all_bn_logn = vcat([bn_samples_by_run[label].logn for label in comparison_run_labels]...)
    all_bn_logB = vcat([bn_samples_by_run[label].logB for label in comparison_run_labels]...)
    bn_nlo, bn_nhi = quantile(all_bn_logn, (0.001, 0.999))
    bn_Blo, bn_Bhi = quantile(all_bn_logB, (0.001, 0.999))
    if !(bn_nhi > bn_nlo)
        bn_nlo, bn_nhi = bn_nlo - 0.5, bn_nhi + 0.5
    end
    if !(bn_Bhi > bn_Blo)
        bn_Blo, bn_Bhi = bn_Blo - 0.5, bn_Bhi + 0.5
    end
    bn_logn_edges = collect(range(bn_nlo, bn_nhi; length = Int(nbins) + 1))
    bn_logB_edges = collect(range(bn_Blo, bn_Bhi; length = Int(nbins) + 1))
    bn_profiles = Dict(label => binned_magnetic_density(
            bn_samples_by_run[label], bn_logn_edges)
        for label in comparison_run_labels)
    bn_fits = Dict(label => magnetic_density_fit(
            bn_profiles[label], Float64(bn_fit_min_density))
        for label in comparison_run_labels)

    active_bn_samples = haskey(bn_samples_by_run, selected_run) ?
        bn_samples_by_run[selected_run] : magnetic_density_samples(cube)
    active_bn_histogram = fit(Histogram,
        (clamp.(active_bn_samples.logn, bn_nlo, prevfloat(bn_nhi)),
            clamp.(active_bn_samples.logB, bn_Blo, prevfloat(bn_Bhi))),
        (bn_logn_edges, bn_logB_edges))
    bn_histogram_total = sum(active_bn_histogram.weights)
    bn_log_probability = map(active_bn_histogram.weights) do weight
        weight > 0 && bn_histogram_total > 0 ? log10(weight / bn_histogram_total) : NaN
    end
    nothing
end

# ╔═╡ b1000003-6f8c-4d0c-9a10-000000000003
begin
    fig_bn = Figure(size = (1100, 470))
    bn_panel = fig_bn[1, 1] = GridLayout()
    bn_joint_axis = latex_axis(bn_panel[1, 1],
        xlabel = L"n\;[\mathrm{cm}^{-3}]", ylabel = L"|B|\;[\mu\mathrm{G}]",
        xscale = log10, yscale = log10,
        xticks = DECADE_TICKS, yticks = DECADE_TICKS,
        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
        xminorticksvisible = true, yminorticksvisible = true)
    bn_heatmap = heatmap!(bn_joint_axis, 10.0 .^ bn_logn_edges,
        10.0 .^ bn_logB_edges, bn_log_probability; colormap = :magma)
    latex_colorbar(bn_panel[1, 2], bn_heatmap,
        label = L"\log_{10}\mathcal{P}_{\mathrm{bin}}")
    colsize!(bn_panel, 2, 22)

    bn_relation_axis = latex_axis(fig_bn[1, 2],
        xlabel = L"n\;[\mathrm{cm}^{-3}]", ylabel = L"|B|\;[\mu\mathrm{G}]",
        xscale = log10, yscale = log10,
        xticks = DECADE_TICKS, yticks = DECADE_TICKS,
        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
        xminorticksvisible = true, yminorticksvisible = true)
    for label in comparison_run_labels
        profile = bn_profiles[label]
        valid = isfinite.(profile.medians)
        x = 10.0 .^ profile.centers[valid]
        median_B = 10.0 .^ profile.medians[valid]
        band!(bn_relation_axis, x, 10.0 .^ profile.lower[valid],
            10.0 .^ profile.upper[valid]; color = (run_colors[label], 0.15))
        lines!(bn_relation_axis, x, median_B; color = run_colors[label], linewidth = 2.8)
        scatter!(bn_relation_axis, x, median_B; color = run_colors[label], markersize = 6)
        fit_result = bn_fits[label]
        if isfinite(fit_result.slope) && !isempty(x)
            fit_x = extrema(x)
            fit_y = 10.0 .^ (fit_result.intercept .+
                fit_result.slope .* log10.(collect(fit_x)))
            lines!(bn_relation_axis, collect(fit_x), fit_y;
                color = run_colors[label], linestyle = :dash, linewidth = 2,
                label = latexstring(legend_run_label(label),
                    raw";\;\kappa=", @sprintf("%.3f", fit_result.slope)))
        end
    end
    axislegend(bn_relation_axis; position = :lt, framevisible = false)
    stable_pluto_figure(display_bn_relation, fig_bn)
end

# ╔═╡ b1000004-6f8c-4d0c-9a10-000000000004
begin
    function hro_products(c, density_bin_count; angle_bin_count = 18)
        local_n = number_density(c.rho)
        log_density = safe_log10.(local_n)
        spacing = c.L ./ size(c.rho)
        gx = periodic_derivative(log_density, 1, spacing[1])
        gy = periodic_derivative(log_density, 2, spacing[2])
        gz = periodic_derivative(log_density, 3, spacing[3])
        gradient_norm = sqrt.(gx .^ 2 .+ gy .^ 2 .+ gz .^ 2)
        field_norm = magnetic_fields(c).B
        valid = isfinite.(log_density) .& isfinite.(gradient_norm) .&
            isfinite.(field_norm) .& (gradient_norm .> 0) .& (field_norm .> 0)
        cosine_gradient = abs.(c.bx .* gx .+ c.by .* gy .+ c.bz .* gz) ./
            max.(field_norm .* gradient_norm, eps(Float64))
        # 0°: B follows an isodensity structure; 90°: B crosses it.
        angles = asind.(clamp.(cosine_gradient[valid], 0.0, 1.0))
        densities = log_density[valid]
        probability_edges = collect(range(0.0, 1.0; length = density_bin_count + 1))
        density_edges = unique(quantile(densities, probability_edges))
        length(density_edges) >= 2 || (density_edges = [minimum(densities), maximum(densities) + eps()])
        angle_edges = collect(range(0.0, 90.0; length = angle_bin_count + 1))
        angle_centers = (angle_edges[1:end-1] .+ angle_edges[2:end]) ./ 2
        bin_count = length(density_edges) - 1
        histograms = zeros(Float64, length(angle_centers), bin_count)
        shape = fill(NaN, bin_count)
        density_centers = fill(NaN, bin_count)
        counts = zeros(Int, bin_count)
        for bin in 1:bin_count
            members = (densities .>= density_edges[bin]) .&
                (bin == bin_count ? densities .<= density_edges[bin + 1] :
                    densities .< density_edges[bin + 1])
            local_angles = angles[members]
            counts[bin] = length(local_angles)
            isempty(local_angles) && continue
            density_centers[bin] = median(densities[members])
            histogram = fit(Histogram, local_angles, angle_edges).weights
            histograms[:, bin] .= histogram ./ max(sum(histogram), 1)
            central = count(angle -> angle < 22.5, local_angles)
            edge = count(angle -> angle > 67.5, local_angles)
            shape[bin] = (central - edge) / max(central + edge, 1)
        end
        (; angle_centers, histograms, density_centers, shape, counts)
    end

    hro_by_run = Dict(label => hro_products(
            comparison_cube(label), Int(hro_density_bin_count))
        for label in comparison_run_labels)
    active_hro = haskey(hro_by_run, selected_run) ? hro_by_run[selected_run] :
        hro_products(cube, Int(hro_density_bin_count))

    fig_hro = Figure(size = (1100, 470))
    hro_hist_axis = latex_axis(fig_hro[1, 1],
        xlabel = L"\phi_{B,\,\mathrm{structure}}\;[{}^\circ]",
        ylabel = L"\mathcal{P}(\phi)")
    representative_bins = unique([1, cld(size(active_hro.histograms, 2), 2),
        size(active_hro.histograms, 2)])
    for (style_index, bin) in enumerate(representative_bins)
        isfinite(active_hro.density_centers[bin]) || continue
        lines!(hro_hist_axis, active_hro.angle_centers,
            active_hro.histograms[:, bin]; color = MHD_COLORS[style_index],
            linewidth = 2.8,
            label = latexstring(raw"n_{\mathrm{med}}=",
                @sprintf("%.3g", 10.0^active_hro.density_centers[bin]),
                raw"\;\mathrm{cm}^{-3}"))
    end
    axislegend(hro_hist_axis; position = :ct, framevisible = false)

    hro_shape_axis = latex_axis(fig_hro[1, 2],
        xlabel = L"n\;[\mathrm{cm}^{-3}]", ylabel = L"\zeta_{\mathrm{HRO}}",
        xscale = log10, xticks = DECADE_TICKS,
        xminorticks = IntervalsBetween(9), xminorticksvisible = true)
    hlines!(hro_shape_axis, [0.0]; color = (:gray45, 0.65), linestyle = :dash)
    for label in comparison_run_labels
        product = hro_by_run[label]
        valid = isfinite.(product.density_centers) .& isfinite.(product.shape)
        lines!(hro_shape_axis, 10.0 .^ product.density_centers[valid],
            product.shape[valid]; color = run_colors[label], linewidth = 2.8,
            label = latexstring(legend_run_label(label)))
        scatter!(hro_shape_axis, 10.0 .^ product.density_centers[valid],
            product.shape[valid]; color = run_colors[label], markersize = 6)
    end
    axislegend(hro_shape_axis; position = :lb, framevisible = false)
    stable_pluto_figure(display_hro, fig_hro)
end

# ╔═╡ b1000005-6f8c-4d0c-9a10-000000000005
begin
    function periodic_gaussian_smooth_2d(image, fwhm_pixels)
        fwhm_pixels <= 0 && return Float64.(image)
        nx, ny = size(image)
        sigma = Float64(fwhm_pixels) / (2sqrt(2log(2)))
        kx = reshape([i <= nx ÷ 2 ? i : i - nx for i in 0:nx-1] ./ nx, nx, 1)
        ky = reshape([j <= ny ÷ 2 ? j : j - ny for j in 0:ny-1] ./ ny, 1, ny)
        transfer = @. exp(-2pi^2 * sigma^2 * (kx^2 + ky^2))
        valid = isfinite.(image)
        filled = ifelse.(valid, Float64.(image), 0.0)
        smoothed = real.(ifft(fft(filled) .* transfer))
        normalization = real.(ifft(fft(Float64.(valid)) .* transfer))
        map((value, weight) -> weight > sqrt(eps(Float64)) ? value / weight : NaN,
            smoothed, normalization)
    end

    function hog_products(c, line_of_sight, plane_dimensions;
            smoothing_fwhm = 2.0, gradient_percentile = 20.0,
            logarithmic_maps = true, angle_bin_count = 18)
        local_n = number_density(c.rho)
        local_B = GAUSS_TO_MICROGAUSS .* magnetic_fields(c).B
        column = finite_sum_dims(local_n, line_of_sight)
        projected_B = weighted_project(local_B, c.rho, line_of_sight)
        image1 = logarithmic_maps ? safe_log10.(column) : Float64.(column)
        image2 = logarithmic_maps ? safe_log10.(projected_B) : Float64.(projected_B)
        image1 = periodic_gaussian_smooth_2d(image1, smoothing_fwhm)
        image2 = periodic_gaussian_smooth_2d(image2, smoothing_fwhm)
        d1x = periodic_derivative(image1, 1, 1.0)
        d1y = periodic_derivative(image1, 2, 1.0)
        d2x = periodic_derivative(image2, 1, 1.0)
        d2y = periodic_derivative(image2, 2, 1.0)
        norm1, norm2 = hypot.(d1x, d1y), hypot.(d2x, d2y)
        positive1 = filter(x -> isfinite(x) && x > 0, vec(norm1))
        positive2 = filter(x -> isfinite(x) && x > 0, vec(norm2))
        threshold1 = isempty(positive1) ? 0.0 : quantile(positive1, gradient_percentile / 100)
        threshold2 = isempty(positive2) ? 0.0 : quantile(positive2, gradient_percentile / 100)
        valid = isfinite.(d1x) .& isfinite.(d1y) .& isfinite.(d2x) .&
            isfinite.(d2y) .& (norm1 .> threshold1) .& (norm2 .> threshold2)
        dot_product = d1x .* d2x .+ d1y .* d2y
        cosine = abs.(dot_product[valid]) ./
            max.(norm1[valid] .* norm2[valid], eps(Float64))
        angles = acosd.(clamp.(cosine, 0.0, 1.0))
        angle_edges = collect(range(0.0, 90.0; length = angle_bin_count + 1))
        angle_centers = (angle_edges[1:end-1] .+ angle_edges[2:end]) ./ 2
        histogram = isempty(angles) ? zeros(length(angle_centers)) :
            Float64.(fit(Histogram, angles, angle_edges).weights)
        sum(histogram) > 0 && (histogram ./= sum(histogram))
        cos2 = cosd.(2 .* angles)
        sin2 = sind.(2 .* angles)
        normalized_prs = isempty(cos2) ? NaN : mean(cos2)
        prs = isempty(cos2) ? NaN : sum(cos2) / sqrt(length(cos2) / 2)
        rvl = isempty(cos2) ? NaN : hypot(mean(cos2), mean(sin2))
        (; angle_centers, histogram, normalized_prs, prs, rvl,
            ngood = length(angles))
    end

    hog_by_run = Dict(label => hog_products(comparison_cube(label), los_dim, sky_dims;
            smoothing_fwhm = Float64(hog_smoothing_fwhm_pix),
            gradient_percentile = Float64(hog_gradient_percentile),
            logarithmic_maps = hog_logarithmic_maps)
        for label in comparison_run_labels)

    fig_hog = Figure(size = (1100, 470))
    hog_hist_axis = latex_axis(fig_hog[1, 1],
        xlabel = L"\phi_{\nabla N_{\mathrm H},\,\nabla |B|}\;[{}^\circ]",
        ylabel = L"\mathcal{P}(\phi)")
    for label in comparison_run_labels
        product = hog_by_run[label]
        lines!(hog_hist_axis, product.angle_centers, product.histogram;
            color = run_colors[label], linewidth = 2.8,
            label = latexstring(legend_run_label(label)))
    end
    axislegend(hog_hist_axis; position = :rt, framevisible = false)

    hog_positions = collect(1:length(comparison_run_labels))
    hog_values = [hog_by_run[label].normalized_prs for label in comparison_run_labels]
    hog_prs_axis = latex_axis(fig_hog[1, 2],
        xlabel = L"\mathrm{simulation}", ylabel = L"V/V_{\max}",
        xticks = hog_positions,
        xtickformat = values -> [latexstring(legend_run_label(
            comparison_run_labels[clamp(round(Int, value), 1,
                length(comparison_run_labels))])) for value in values])
    hlines!(hog_prs_axis, [0.0]; color = (:gray45, 0.65), linestyle = :dash)
    barplot!(hog_prs_axis, hog_positions, hog_values;
        color = [run_colors[label] for label in comparison_run_labels])
    stable_pluto_figure(display_hog, fig_hog)
end

# ╔═╡ 298bd579-bb28-48b7-8c55-ea74804b9837
md"""
---

## 8. Energy partition by density and time

### Density-binned energy ratios

The figure compares three energy ratios as functions of physical number density. Color identifies the simulation run and line style identifies the snapshot. The horizontal reference at unity marks equal energies; values above unity indicate that the numerator dominates.

Each point is a ratio of energies summed within one density bin, rather than a mean of cell-by-cell ratios. The CGS energy densities are $E_{\mathrm{kin}}=\rho|\delta\mathbf v|^2/2$, $E_{\mathrm{mag}}=B^2/(8\pi)$, and $E_{\mathrm{therm}}=P/(\gamma-1)$ in $\mathrm{erg\,cm}^{-3}$. For $\gamma=1$, the notebook adopts the isothermal convention $E_{\mathrm{therm}}=P$. Cubes selected for the comparative profiles are read and reduced to these per-cell quantities in a separate Pluto cell; changing $\gamma$ or the density bins therefore rebins the cached arrays without rereading the RAMSES files.

**Display energy ratios by density:** $(@bind display_energy_ratios PlutoUI.CheckBox(default = true))

**Display energy ratios versus time:** $(@bind display_energy_time PlutoUI.CheckBox(default = false))

| Setting | Control |
|:--|:--|
| Snapshots shown in energy reports | $(@bind energy_snapshot_indices PlutoUI.MultiSelect(collect(1:maximum_snapshot_count); default = [min(Int(selected_snapshot), maximum_snapshot_count)])) |
| $E_{\mathrm{kin}}/E_{\mathrm{mag}}$ | $(@bind show_energy_kin_mag PlutoUI.CheckBox(default = true)) |
| $E_{\mathrm{therm}}/E_{\mathrm{mag}}$ | $(@bind show_energy_therm_mag PlutoUI.CheckBox(default = true)) |
| $E_{\mathrm{kin}}/E_{\mathrm{therm}}$ | $(@bind show_energy_kin_therm PlutoUI.CheckBox(default = true)) |

### Time evolution of the energy ratios

The same ratio checkboxes select the time-evolution panels. Each curve uses energies integrated over the complete simulation volume. Color identifies the run, and the horizontal reference at unity marks equipartition.
"""

# ╔═╡ e4a64ef5-9e6e-4ae2-8daf-42e4fc1724ba
begin
    function energy_density_samples(c)
        m = magnetic_fields(c)
        u = turbulent_velocity(c)
        (
            log_density = vec(safe_log10.(number_density(c.rho))),
            emag = vec(m.B2 ./ (8pi)),
            ekin = vec(0.5 .* c.rho .* u.dv2 .* KM_CM^2),
            pressure = vec(Float64.(c.P)),
        )
    end

    function energy_ratios_by_density(samples, edges, gamma)
        s = samples.log_density
        Emag = samples.emag
        Ekin = samples.ekin
        pressure = samples.pressure
        thermal_factor = gamma > 1 + sqrt(eps(Float64)) ? 1 / (gamma - 1) : 1.0
        idx = clamp.(searchsortedlast.(Ref(edges), s), 1, length(edges) - 1)
        emag = zeros(length(edges) - 1)
        ekin = zeros(length(edges) - 1)
        etherm = zeros(length(edges) - 1)
        for i in eachindex(idx)
            b = idx[i]
            Etherm = thermal_factor * pressure[i]
            if isfinite(s[i]) && isfinite(Emag[i]) && isfinite(Ekin[i]) &&
                    isfinite(Etherm) && edges[b] <= s[i] <= edges[b + 1]
                emag[b] += Emag[i]
                ekin[b] += Ekin[i]
                etherm[b] += Etherm
            end
        end
        positive_ratio(numerator, denominator) = [
            numerator[i] > 0 && denominator[i] > 0 ? numerator[i] / denominator[i] : NaN
            for i in eachindex(numerator)
        ]
        (
            kin_mag = positive_ratio(ekin, emag),
            therm_mag = positive_ratio(etherm, emag),
            kin_therm = positive_ratio(ekin, etherm),
            mag_kin = positive_ratio(emag, ekin),
        )
    end

end

# ╔═╡ a1bf5e62-ecb8-47d7-8692-c33929c831ac
selected_energy_samples = energy_density_samples(cube)

# ╔═╡ 5983ddd0-3ea9-4304-a1e2-390065309aef
begin
    density_values = finite_values(logn)
    isempty(density_values) && error("No finite positive density is available.")
    density_lo, density_hi = quantile(density_values, (0.001, 0.999))
    density_lo == density_hi && ((density_lo, density_hi) =
        (density_lo - 0.5, density_hi + 0.5))
    density_edges = range(density_lo, density_hi; length = nbins + 1)
    density_centers = (density_edges[1:end-1] .+ density_edges[2:end]) ./ 2
    density_number_centers = 10.0 .^ density_centers
    selected_energy_ratios = energy_ratios_by_density(
        selected_energy_samples, density_edges, gamma)
end

# ╔═╡ b7c04a6f-b24d-48ca-97d2-5a77d297aa8d
begin
    valid_energy_snapshots = sort(unique(filter(i -> i >= 1, Int.(energy_snapshot_indices))))
    energy_profiles = Dict{Tuple{String, Int}, Any}()
    for label in comparison_run_labels
        for snapshot_index in valid_energy_snapshots
            snapshot_index <= length(run_files[label]) || continue
            # Build and retain only the one-dimensional profile. The cube from
            # this iteration becomes collectible before the next file is read.
            local_cube = label == selected_run && snapshot_index == selected_snapshot ?
                cube : load_cube(run_files[label][snapshot_index])
            local_samples = energy_density_samples(local_cube)
            energy_profiles[(label, snapshot_index)] = energy_ratios_by_density(
                local_samples, density_edges, gamma)
        end
    end
end

# ╔═╡ 2c4762f5-4780-4bfa-a320-e4102d9f5d6c
nothing

# ╔═╡ 72b2e359-c054-4efc-a834-a2f5fa99e3fc
begin
    energy_panel_keys = Symbol[]
    show_energy_kin_mag && push!(energy_panel_keys, :kin_mag)
    show_energy_therm_mag && push!(energy_panel_keys, :therm_mag)
    show_energy_kin_therm && push!(energy_panel_keys, :kin_therm)
    energy_line_styles = [:solid, :dash, :dot, :dashdot]
    snapshot_styles = Dict(
        snapshot_index => energy_line_styles[mod1(style_index, length(energy_line_styles))]
        for (style_index, snapshot_index) in enumerate(valid_energy_snapshots)
    )
    if isempty(energy_panel_keys)
        fig_energy = Figure(size = (900, 180))
        Label(fig_energy[1, 1], L"\mathrm{Select\ at\ least\ one\ energy\ ratio.}", fontsize = 20)
    else
        fig_energy = Figure(size = (500length(energy_panel_keys), 520))
        energy_axes = Dict{Symbol, Any}()
        for (index, key) in enumerate(energy_panel_keys)
            ylabel = key == :kin_mag ? L"E_{\mathrm{kin}}/E_{\mathrm{mag}}" :
                key == :therm_mag ? L"E_{\mathrm{therm}}/E_{\mathrm{mag}}" :
                L"E_{\mathrm{kin}}/E_{\mathrm{therm}}"
            energy_axes[key] = latex_axis(fig_energy[1, index],
                xlabel = L"n\;[\mathrm{cm}^{-3}]", ylabel = ylabel,
                xscale = log10, yscale = log10,
                xticks = DECADE_TICKS, yticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                xminorticksvisible = true, yminorticksvisible = true)
            hlines!(energy_axes[key], [1.0]; color = (:gray45, 0.55),
                linestyle = :dash, linewidth = 1.5)
        end
        for label in comparison_run_labels, snapshot_index in valid_energy_snapshots
            haskey(energy_profiles, (label, snapshot_index)) || continue
            profile = energy_profiles[(label, snapshot_index)]
            style = snapshot_styles[snapshot_index]
            for key in energy_panel_keys
                lines!(energy_axes[key], density_number_centers, getfield(profile, key);
                    color = run_colors[label], linestyle = style, linewidth = 2.5)
            end
        end
        run_legend_elements = [LineElement(color = run_colors[label], linewidth = 2.5) for label in comparison_run_labels]
        snapshot_legend_elements = [
            LineElement(color = :black, linestyle = snapshot_styles[index], linewidth = 2.5)
            for index in valid_energy_snapshots
        ]
        energy_run_labels = legend_run_label.(comparison_run_labels)
        energy_snapshot_labels = ["Snapshot $(index)" for index in valid_energy_snapshots]
        Legend(fig_energy[2, 1:length(energy_panel_keys)],
            [run_legend_elements, snapshot_legend_elements],
            [energy_run_labels, energy_snapshot_labels],
            ["Run", "Snapshot"];
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    display_energy_ratios ? fig_energy : nothing
end

# ╔═╡ 4e8ac6d4-4ee4-4e41-a63a-1c93213da33f
begin
    energy_time_specs = NamedTuple[]
    show_energy_kin_mag && push!(energy_time_specs,
        (field = :kin_mag, ylabel = L"E_{\mathrm{kin}}/E_{\mathrm{mag}}"))
    show_energy_therm_mag && push!(energy_time_specs,
        (field = :therm_mag, ylabel = L"E_{\mathrm{therm}}/E_{\mathrm{mag}}"))
    show_energy_kin_therm && push!(energy_time_specs,
        (field = :kin_therm, ylabel = L"E_{\mathrm{kin}}/E_{\mathrm{therm}}"))

    if isempty(energy_time_specs)
        fig_energy_time = Figure(size = (900, 180))
        Label(fig_energy_time[1, 1],
            L"\mathrm{Select\ at\ least\ one\ time\ dependent\ energy\ ratio.}", fontsize = 20)
    else
        fig_energy_time = Figure(size = (500length(energy_time_specs), 470))
        for (index, spec) in enumerate(energy_time_specs)
            ax = latex_axis(fig_energy_time[1, index],
                xlabel = L"t\;[\mathrm{Myr}]", ylabel = spec.ylabel, yscale = log10,
                yticks = DECADE_TICKS,
                yminorticks = IntervalsBetween(9), yminorticksvisible = true)
            hlines!(ax, [1.0]; color = (:gray45, 0.55), linestyle = :dash, linewidth = 1.5)
            for label in comparison_run_labels
                series = all_series[label]
                times = Float64.(getfield.(series, :t))
                ratio = Float64.(getfield.(series, spec.field))
                lines!(ax, times, ratio; color = run_colors[label], linewidth = 2.6)
                scatter!(ax, times, ratio; color = run_colors[label], markersize = 6)
            end
        end
        energy_time_legend_elements = [
            LineElement(color = run_colors[label], linewidth = 2.6) for label in comparison_run_labels
        ]
        energy_time_legend_labels = legend_run_label.(comparison_run_labels)
        Legend(fig_energy_time[2, 1:length(energy_time_specs)],
            energy_time_legend_elements, energy_time_legend_labels;
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    display_energy_time ? fig_energy_time : nothing
end

# ╔═╡ d87379b7-3527-45a8-bc60-bec191c499af
md"""
---

## 9. Vorticity and enstrophy

The first map shows the line-of-sight mean of the vorticity magnitude, $\langle|\boldsymbol\omega|\rangle_{\mathrm{LOS}}$. The second shows the line-of-sight mean enstrophy,

```math
\mathcal{E}_{\omega}=\frac{1}{2}|\boldsymbol\omega|^2,
```

in $\mathrm{Myr}^{-2}$. Both maps use the active run, snapshot, and line of sight.

**Display vorticity and enstrophy maps:** $(@bind display_vorticity_figure PlutoUI.CheckBox(default = true))  
**Display density-binned enstrophy diagnostics:** $(@bind display_enstrophy_density PlutoUI.CheckBox(default = true))

| Vorticity and enstrophy diagnostic | Display |
|:--|:--:|
| Vorticity heatmap | $(@bind show_vorticity_map PlutoUI.CheckBox(default = true)) |
| Enstrophy heatmap | $(@bind show_enstrophy_map PlutoUI.CheckBox(default = true)) |
| Enstrophy profiles by density bin | $(@bind show_enstrophy_density_profiles PlutoUI.CheckBox(default = true)) |
| Density-bin weighting | $(@bind enstrophy_density_weighting PlutoUI.Select(["volume", "mass"]; default = "volume")) |
"""

# ╔═╡ 3190e127-1d53-49f1-bfab-b9645910c2c6
begin
    function isotropic_power_spectrum(components, box_length_pc; prefactor = 1.0)
        Fpower = zeros(Float64, size(first(components)))
        for component in components
            replacement = finite_mean(component; default = 0.0)
            finite_component = ifelse.(isfinite.(component), Float64.(component), replacement)
            Fpower .+= abs2.(fft(finite_component))
        end
        n = size(Fpower)
        ntot = length(Fpower)
        kmax = floor(Int, sqrt(sum((ni ÷ 2)^2 for ni in n)))
        shell_power = zeros(kmax + 1)
        shell_count = zeros(Int, kmax + 1)
        fft_mode(i, ni) = (i - 1 <= ni ÷ 2) ? i - 1 : i - 1 - ni
        for I in CartesianIndices(Fpower)
            kmag = sqrt(sum(fft_mode(I[d], n[d])^2 for d in 1:3))
            shell = round(Int, kmag) + 1
            shell <= length(shell_power) || continue
            shell_power[shell] += prefactor * Fpower[I] / ntot^2
            shell_count[shell] += 1
        end
        modes = collect(0:kmax)
        valid = (modes .> 0) .& (shell_count .> 0) .& (shell_power .> 0)
        dk = 2pi / box_length_pc
        modes[valid] .* dk, shell_power[valid] ./ dk
    end

    spectrum_box_length_pc = cbrt(prod(cube.L))
    spectrum_dk = 2pi / spectrum_box_length_pc
    spectrum_low_k_limit = 3spectrum_dk
    krho, Prho = isotropic_power_spectrum((number_density_cells,), spectrum_box_length_pc)
    kv, Pv = isotropic_power_spectrum((turb.dvx, turb.dvy, turb.dvz), spectrum_box_length_pc; prefactor = 0.5)
    komega, Pomega = isotropic_power_spectrum((omega.wx, omega.wy, omega.wz), spectrum_box_length_pc)
    # The shell power is quadratic in the field, so squaring the Gauss-to-microgauss
    # factor into the prefactor puts E_B in microgauss^2 pc without scaling the
    # three components into new cube-sized arrays first.
    kb, Pb = isotropic_power_spectrum((cube.bx, cube.by, cube.bz), spectrum_box_length_pc;
        prefactor = 0.5GAUSS_TO_MICROGAUSS^2)
    spectrum_maximum_k = maximum(vcat(krho, kv, komega, kb))
    spectrum_k_choices = spectrum_dk:spectrum_dk:spectrum_maximum_k
    spectrum_default_k_min = min(spectrum_low_k_limit, spectrum_maximum_k)
    spectrum_default_k_max = max(spectrum_default_k_min,
        spectrum_dk * floor(Int, spectrum_maximum_k / (2spectrum_dk)))
end

# ╔═╡ df880704-b12e-49eb-84f2-67b6f9583a8a
begin
    vorticity_map_specs = NamedTuple[]
    show_vorticity_map && push!(vorticity_map_specs,
        (data = safe_log10.(omega_map), colormap = :inferno,
            label = L"\log_{10}\!\left(\langle|\omega|\rangle/\mathrm{Myr}^{-1}\right)"))
    show_enstrophy_map && push!(vorticity_map_specs,
        (data = safe_log10.(enstrophy_map), colormap = :magma,
            label = L"\log_{10}\!\left(\langle\mathcal{E}_{\omega}\rangle/\mathrm{Myr}^{-2}\right)"))
    if isempty(vorticity_map_specs)
        fig_vorticity = Figure(size = (900, 180))
        Label(fig_vorticity[1, 1],
            L"\mathrm{Select\ at\ least\ one\ vorticity\ or\ enstrophy\ map.}", fontsize = 20)
    else
        fig_vorticity = Figure(size = (620length(vorticity_map_specs), 520))
        for (index, spec) in enumerate(vorticity_map_specs)
            panel = fig_vorticity[1, index] = GridLayout()
            axis = latex_axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            heat = heatmap!(axis, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap)
            latex_colorbar(panel[1, 2], heat; label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_vorticity_figure ? fig_vorticity : nothing
end

# ╔═╡ b02ceda6-0a0d-4543-91f8-a6669231ec72
begin
    function enstrophy_by_density(c, edges, weighting)
        local_omega = vorticity(c)
        local_enstrophy = 0.5 .* local_omega.magnitude .^ 2
        log_density = vec(safe_log10.(number_density(c.rho)))
        enstrophy_values = vec(local_enstrophy)
        weights = weighting == "mass" ? vec(Float64.(c.rho)) : ones(length(c.rho))
        bin_index = clamp.(searchsortedlast.(Ref(edges), log_density), 1, length(edges) - 1)
        weighted_sum = zeros(Float64, length(edges) - 1)
        weight_sum = zeros(Float64, length(edges) - 1)
        for cell in eachindex(bin_index)
            bin = bin_index[cell]
            if edges[bin] <= log_density[cell] <= edges[bin + 1] &&
                    isfinite(enstrophy_values[cell]) && isfinite(weights[cell])
                weighted_sum[bin] += weights[cell] * enstrophy_values[cell]
                weight_sum[bin] += weights[cell]
            end
        end
        [weight_sum[index] > 0 ? weighted_sum[index] / weight_sum[index] : NaN
            for index in eachindex(weight_sum)]
    end

    enstrophy_snapshot_indices = Dict(
        label => min(Int(selected_snapshot), length(run_files[label])) for label in comparison_run_labels)
    enstrophy_profiles = Dict{String, Any}()
    enstrophy_mach_by_run = Dict{String, Float64}()
    for label in comparison_run_labels
        local_cube = comparison_cube(label)
        enstrophy_profiles[label] = enstrophy_by_density(
            local_cube, density_edges, enstrophy_density_weighting)
        comparison_kind == :mach && (enstrophy_mach_by_run[label] =
            bulk_metrics_from_cube(local_cube, Float64(gamma)).mach)
    end
    function enstrophy_parameter_label(label)
        if comparison_kind == :mach
            mach_value = enstrophy_mach_by_run[label]
            return latexstring("\\mathcal{M}=",
                @sprintf("%.3g", mach_value))
        end
        latex_run_label(label)
    end
end

# ╔═╡ 873f7ef2-719b-4ae6-b015-1a23c6c27836
begin
    if !show_enstrophy_density_profiles
        fig_enstrophy_density = Figure(size = (900, 180))
        Label(fig_enstrophy_density[1, 1],
            L"\mathrm{Density-binned\ enstrophy\ profiles\ are\ disabled.}",
            fontsize = 20)
    else
        fig_enstrophy_density = Figure(size = (680, 500))
        axis = latex_axis(fig_enstrophy_density[1, 1],
            xlabel = L"n\;[\mathrm{cm}^{-3}]",
            ylabel = L"\langle\mathcal{E}_{\omega}\rangle_n\;[\mathrm{Myr}^{-2}]",
            xscale = log10, yscale = log10,
            xticks = DECADE_TICKS, yticks = DECADE_TICKS,
            xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
            xminorticksvisible = true, yminorticksvisible = true)
        for label in comparison_run_labels
            profile = enstrophy_profiles[label]
            valid = isfinite.(profile) .& (profile .> 0)
            lines!(axis, density_number_centers[valid], profile[valid];
                color = run_colors[label], linewidth = 2.5,
                label = enstrophy_parameter_label(label))
            scatter!(axis, density_number_centers[valid], profile[valid];
                color = run_colors[label], markersize = 5)
        end
        axislegend(axis; position = :lt, framevisible = false)
    end
    display_enstrophy_density ? fig_enstrophy_density : nothing
end

# ╔═╡ a8558c31-7dcf-433e-9950-a59e9acf158b
md"""
---

## 10. Isotropic power spectra

The panels show number-density, turbulent-velocity, vorticity, and magnetic spectra on base-10 logarithmic axes in both $k$ and spectral power. Fourier power is summed in spherical shells and normalized consistently with Parseval's theorem. Dividing shell power by $\Delta k=2\pi/L$ gives spectral densities satisfying $\int E_f(k)\,\mathrm{d}k=\langle|f|^2\rangle$, with the stated velocity prefactor. Velocity and vorticity include all three vector components, and $k$ is expressed in $\mathrm{pc}^{-1}$.

The first shells contain few Fourier modes and are correspondingly noisy. The shaded region $k<3\,\Delta k$ marks these box-scale modes; it should be excluded when estimating an inertial-range slope. The threshold is a visual reliability guide, not a claim that an inertial range necessarily begins at $3\,\Delta k$.

The fitted model is $E(k)=A k^\alpha$ over the selected interval. Each optional reference slope is normalized to the fitted spectrum at the geometric center of that interval.

The magnetic panel shows $E_B(k)=\tfrac12\langle|\mathbf B|^2\rangle_k$ in $\mu\mathrm G^2\,\mathrm{pc}$, so it is an energy spectrum up to the constant $1/4\pi$.

Two different references are offered, because they describe different fields:

- **Kolmogorov**, $\alpha=-5/3$, for the density, velocity, and vorticity panels;
- **Kazantsev** (1968), $\alpha=+3/2$, for the magnetic panel.

The Kazantsev slope is the prediction of the *kinematic* small-scale dynamo: while the field is still too weak to react back on the flow, magnetic energy piles up at small scales and $E_B(k)\propto k^{3/2}$ between the forcing scale and the resistive scale. Its **positive** exponent is the signature of that regime — magnetic energy peaks at the *smallest* resolved scales, not the largest.

This reference is therefore only meaningful while the dynamo is still kinematic. Once the field saturates, back-reaction flattens and then bends the spectrum, and a $k^{3/2}$ fit stops being informative. Use the exponential-growth window fitted in section 5 to decide whether the selected snapshot is still in the kinematic phase.

**Display density, velocity, vorticity, and magnetic power spectra:** $(@bind display_power_spectra PlutoUI.CheckBox(default = true))

| Power-spectrum panel | Display |
|:--|:--:|
| Number-density spectrum | $(@bind show_spectrum_density PlutoUI.CheckBox(default = true)) |
| Velocity spectrum | $(@bind show_spectrum_velocity PlutoUI.CheckBox(default = true)) |
| Vorticity spectrum | $(@bind show_spectrum_vorticity PlutoUI.CheckBox(default = true)) |
| Magnetic spectrum | $(@bind show_spectrum_magnetic PlutoUI.CheckBox(default = true)) |
| Display fitted slopes | $(@bind show_spectrum_slopes PlutoUI.CheckBox(default = true)) |
| Minimum fitted wavenumber [$\mathrm{pc}^{-1}$] | $(@bind spectrum_fit_k_min PlutoUI.NumberField(spectrum_k_choices; default = spectrum_default_k_min)) |
| Maximum fitted wavenumber [$\mathrm{pc}^{-1}$] | $(@bind spectrum_fit_k_max PlutoUI.NumberField(spectrum_k_choices; default = spectrum_default_k_max)) |
| Display the Kolmogorov $k^{-5/3}$ reference | $(@bind show_kolmogorov_spectrum PlutoUI.CheckBox(default = true)) |
| Display the Kazantsev $k^{3/2}$ reference | $(@bind show_kazantsev_spectrum PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 3a731972-3404-478c-a572-00a05ab652b1
begin
    function power_law_slope(k, power, minimum_k, maximum_k)
        lower, upper = minmax(Float64(minimum_k), Float64(maximum_k))
        valid = isfinite.(k) .& isfinite.(power) .& (k .> 0) .& (power .> 0) .&
            (k .>= lower) .& (k .<= upper)
        x, y = log10.(Float64.(k[valid])), log10.(Float64.(power[valid]))
        length(x) >= 2 || return (slope = NaN, intercept = NaN, r2 = NaN,
            count = length(x), lower = lower, upper = upper)
        xmean, ymean = mean(x), mean(y)
        denominator = sum(abs2, x .- xmean)
        denominator > 0 || return (slope = NaN, intercept = NaN, r2 = NaN,
            count = length(x), lower = lower, upper = upper)
        slope = sum((x .- xmean) .* (y .- ymean)) / denominator
        intercept = ymean - slope * xmean
        prediction = intercept .+ slope .* x
        residual = sum(abs2, y .- prediction)
        total = sum(abs2, y .- ymean)
        r2 = total > 0 ? 1 - residual / total : NaN
        (; slope, intercept, r2, count = length(x), lower, upper)
    end

    # Turbulent cascade panels quote Kolmogorov; the magnetic panel quotes
    # Kazantsev, whose exponent is positive because kinematic small-scale dynamo
    # action piles magnetic energy up at the smallest scales.
    kolmogorov_reference = show_kolmogorov_spectrum ?
        (exponent = -5 / 3, label = L"k^{-5/3}") : nothing
    kazantsev_reference = show_kazantsev_spectrum ?
        (exponent = 3 / 2, label = L"k^{3/2}") : nothing

    spectrum_specs = NamedTuple[]
    show_spectrum_density && push!(spectrum_specs,
        (k = krho, power = Prho, ylabel = L"E_n(k)\;[\mathrm{cm}^{-6}\,\mathrm{pc}]",
            color = MHD_COLORS[1], reference = kolmogorov_reference))
    show_spectrum_velocity && push!(spectrum_specs,
        (k = kv, power = Pv,
            ylabel = L"E_v(k)\;[(\mathrm{km\,s}^{-1})^2\,\mathrm{pc}]",
            color = MHD_COLORS[3], reference = kolmogorov_reference))
    show_spectrum_vorticity && push!(spectrum_specs,
        (k = komega, power = Pomega,
            ylabel = L"E_\omega(k)\;[\mathrm{Myr}^{-2}\,\mathrm{pc}]",
            color = MHD_COLORS[5], reference = kolmogorov_reference))
    show_spectrum_magnetic && push!(spectrum_specs,
        (k = kb, power = Pb,
            ylabel = L"E_B(k)\;[\mu\mathrm{G}^2\,\mathrm{pc}]",
            color = MHD_COLORS[4], reference = kazantsev_reference))
    if isempty(spectrum_specs)
        fig_spectra = Figure(size = (900, 180))
        Label(fig_spectra[1, 1], L"\mathrm{Select\ at\ least\ one\ power\ spectrum.}", fontsize = 20)
    else
        fig_spectra = Figure(size = (400length(spectrum_specs), 390))
        for (index, spec) in enumerate(spectrum_specs)
            axis = latex_axis(fig_spectra[1, index], xlabel = L"k\;[\mathrm{pc}^{-1}]",
                ylabel = spec.ylabel, xscale = log10, yscale = log10,
                xticks = DECADE_TICKS, yticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                xminorticksvisible = true, yminorticksvisible = true,
                xminorticksize = 4, yminorticksize = 4)
            vspan!(axis, spectrum_dk, spectrum_low_k_limit;
                color = (:gray55, 0.18))
            vlines!(axis, [spectrum_low_k_limit];
                color = (:gray35, 0.75), linestyle = :dash, linewidth = 1.3)
            valid = isfinite.(spec.k) .& isfinite.(spec.power) .&
                (spec.k .> 0) .& (spec.power .> 0)
            lines!(axis, spec.k[valid], spec.power[valid]; color = spec.color,
                linewidth = 2.5, label = L"E(k)")
            scatter!(axis, spec.k[valid], spec.power[valid]; color = spec.color, markersize = 5)
            fit_result = power_law_slope(spec.k, spec.power,
                spectrum_fit_k_min, spectrum_fit_k_max)
            if isfinite(fit_result.slope)
                fit_k = 10.0 .^ range(log10(fit_result.lower),
                    log10(fit_result.upper); length = 100)
                fit_power = 10.0 .^ (fit_result.intercept .+
                    fit_result.slope .* log10.(fit_k))
                if show_spectrum_slopes
                    vspan!(axis, fit_result.lower, fit_result.upper;
                        color = (MHD_COLORS[2], 0.08))
                    lines!(axis, fit_k, fit_power; color = MHD_COLORS[2],
                        linewidth = 2.5, linestyle = :dash,
                        label = latexstring(raw"\alpha=", @sprintf("%.3f", fit_result.slope),
                            raw",\;R^2=", @sprintf("%.3f", fit_result.r2)))
                end
                if !isnothing(spec.reference)
                    pivot_k = sqrt(fit_result.lower * fit_result.upper)
                    pivot_power = 10.0^(fit_result.intercept +
                        fit_result.slope * log10(pivot_k))
                    reference_power = pivot_power .*
                        (fit_k ./ pivot_k) .^ spec.reference.exponent
                    lines!(axis, fit_k, reference_power; color = :black,
                        linewidth = 2.2, linestyle = :dot,
                        label = spec.reference.label)
                end
                (show_spectrum_slopes || !isnothing(spec.reference)) &&
                    axislegend(axis; position = :lb, framevisible = false)
            end
        end
    end
    display_power_spectra ? fig_spectra : nothing
end

# ╔═╡ 24e60849-1c70-4df3-bd17-57d29949b7a6
md"""
---

## 11. Real-space structure functions

The panels measure velocity, vorticity, and magnetic-field increments as functions of spatial separation. For a vector field $\mathbf f$, $S_p^f(\ell)=\left\langle\left|\mathbf f(\mathbf x+\boldsymbol\ell)-\mathbf f(\mathbf x)\right|^p\right\rangle$.

The average includes every cell and shifts along the three periodic grid axes, then averages those three results. This is an axis-sampled estimator of the full vector-increment magnitude: it mixes longitudinal and transverse contributions and is **not** an angularly isotropic average. Longitudinal or transverse intermittency exponents should therefore be measured with dedicated projected increments rather than inferred from these panels. Separation is measured in parsecs; vertical units are the corresponding physical field units raised to order $p$. **Order $p$** selects the increment moment, while **Number of separation samples** controls the balance between radial resolution and computation time.

The FFT spectra and structure-function arrays live in calculation cells separate from their plotting cells, so display checkboxes and other cosmetic plot controls do not rerun the heavy transforms. Changing the active cube, the structure-function order, or the separation sampling does require a new calculation. On grids of $256^3$ cells or larger this can still take appreciable time and memory.

**Display velocity, vorticity, and magnetic structure functions:** $(@bind display_structure_functions PlutoUI.CheckBox(default = true))

| Setting | Control |
|:--|:--|
| Order $p$ | $(@bind structure_order PlutoUI.Slider(1:4; default = 2, show_value = true)) |
| Number of separation samples | $(@bind structure_samples PlutoUI.Slider(6:2:24; default = 12, show_value = true)) |
| Velocity structure function | $(@bind show_structure_velocity PlutoUI.CheckBox(default = true)) |
| Vorticity structure function | $(@bind show_structure_vorticity PlutoUI.CheckBox(default = true)) |
| Magnetic structure function | $(@bind show_structure_magnetic PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 4b16d83f-4a1d-49e7-9270-8f573fd46835
begin
    function vector_structure_function(components, lags, order)
        values = zeros(Float64, length(lags))
        shifted = zeros(Float64, size(first(components)))
        increment2 = similar(shifted)
        for (lag_index, lag) in pairs(lags)
            directional_sum = 0.0
            for dimension in 1:3
                shift = ntuple(d -> d == dimension ? lag : 0, 3)
                fill!(increment2, 0.0)
                for component in components
                    circshift!(shifted, component, shift)
                    @. increment2 += abs2(shifted - component)
                end
                finite_sum = 0.0
                finite_count = 0
                for squared_increment in increment2
                    moment = squared_increment^(order / 2)
                    if isfinite(moment)
                        finite_sum += moment
                        finite_count += 1
                    end
                end
                directional_sum += finite_count > 0 ? finite_sum / finite_count : 0.0
            end
            values[lag_index] = directional_sum / 3
        end
        values
    end

    maximum_lag = max(1, minimum(size(cube.rho)) ÷ 2)
    structure_lags = unique(round.(Int, exp.(range(log(1.0), log(Float64(maximum_lag)); length = structure_samples))))
    structure_separations_pc = structure_lags .* minimum(cube.L ./ size(cube.rho))
    Sv = vector_structure_function((turb.dvx, turb.dvy, turb.dvz), structure_lags, structure_order)
    Somega = vector_structure_function((omega.wx, omega.wy, omega.wz), structure_lags, structure_order)
    SB = vector_structure_function(
        (GAUSS_TO_MICROGAUSS .* cube.bx, GAUSS_TO_MICROGAUSS .* cube.by, GAUSS_TO_MICROGAUSS .* cube.bz),
        structure_lags, structure_order)
end

# ╔═╡ d69dd1ce-312a-48e0-a478-b470b299ed1b
begin
    structure_specs = NamedTuple[]
    show_structure_velocity && push!(structure_specs,
        (values = Sv,
            ylabel = latexstring("S_{", structure_order,
                "}^{v}(\\ell)\\;[(\\mathrm{km\\,s}^{-1})^{", structure_order, "}]"),
            color = MHD_COLORS[3]))
    show_structure_vorticity && push!(structure_specs,
        (values = Somega,
            ylabel = latexstring("S_{", structure_order,
                "}^{\\omega}(\\ell)\\;[\\mathrm{Myr}^{-", structure_order, "}]"),
            color = MHD_COLORS[5]))
    show_structure_magnetic && push!(structure_specs,
        (values = SB,
            ylabel = latexstring("S_{", structure_order,
                "}^{B}(\\ell)\\;[(\\mu\\mathrm{G})^{", structure_order, "}]"),
            color = MHD_COLORS[4]))
    if isempty(structure_specs)
        fig_structure = Figure(size = (900, 180))
        Label(fig_structure[1, 1], L"\mathrm{Select\ at\ least\ one\ structure\ function.}", fontsize = 20)
    else
        fig_structure = Figure(size = (400length(structure_specs), 390))
        for (index, spec) in enumerate(structure_specs)
            axis = latex_axis(fig_structure[1, index], xlabel = L"\ell\;[\mathrm{pc}]",
                ylabel = spec.ylabel, xscale = log10, yscale = log10,
                xticks = DECADE_TICKS, yticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                xminorticksvisible = true, yminorticksvisible = true)
            lines!(axis, structure_separations_pc, spec.values;
                color = spec.color, linewidth = 2.5)
            scatter!(axis, structure_separations_pc, spec.values;
                color = spec.color, markersize = 6)
            x_limits = enclosing_decade_limits(structure_separations_pc)
            y_limits = enclosing_decade_limits(spec.values)
            isnothing(x_limits) || xlims!(axis, x_limits...)
            isnothing(y_limits) || ylims!(axis, y_limits...)
        end
    end
    display_structure_functions ? fig_structure : nothing
end

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

# ╔═╡ d6a2f4b1-59ac-4e77-a10a-4b74c0d89231
md"""
---

## 13. Thermal-dust polarization

Optically thin thermal-dust Stokes emission. Statistical plots compare all selected simulations with shared bins.

| Dust figure | Display |
|:--|:--:|
| Polarization maps | $(@bind display_dust_maps PlutoUI.CheckBox(default = true)) |
| Polarization statistics | $(@bind display_dust_statistics PlutoUI.CheckBox(default = true)) |
| Pixel $I_\nu$, $Q_\nu$, and $U_\nu$ spectra | $(@bind display_dust_pixel_spectrum PlutoUI.CheckBox(default = true)) |

| Dust setting | Control |
|:--|:--|
| Observing frequency [$\mathrm{GHz}$] | $(@bind dust_frequency_GHz PlutoUI.NumberField(30.0:1.0:1200.0; default = 353.0)) |
| Dust temperature [$\mathrm{K}$] | $(@bind dust_temperature_K PlutoUI.NumberField(2.8:0.1:100.0; default = 19.6)) |
| Cross-section at $353\,\mathrm{GHz}$ [$\mathrm{cm^2\,H^{-1}}$] | $(@bind dust_sigma353_cm2 PlutoUI.NumberField(default = 1.0e-26)) |
| Emissivity index $\beta_{\mathrm d}$ | $(@bind dust_beta PlutoUI.NumberField(0.0:0.05:4.0; default = 1.6)) |
| Intrinsic polarization fraction $p_0$ | $(@bind dust_p0 PlutoUI.NumberField(0.0:0.005:0.5; default = 0.20)) |
| Gas mass per H nucleon [$m_H$] | $(@bind dust_mu_H PlutoUI.NumberField(1.0:0.01:2.0; default = 1.4)) |
| $I_\nu$ map | $(@bind show_dust_I PlutoUI.CheckBox(default = true)) |
| $Q_\nu$ map | $(@bind show_dust_Q PlutoUI.CheckBox(default = false)) |
| $U_\nu$ map | $(@bind show_dust_U PlutoUI.CheckBox(default = false)) |
| Polarized intensity $P_\nu$ | $(@bind show_dust_P PlutoUI.CheckBox(default = true)) |
| Polarization fraction $p$ | $(@bind show_dust_fraction PlutoUI.CheckBox(default = true)) |
| Polarization angle $\psi$ | $(@bind show_dust_angle PlutoUI.CheckBox(default = false)) |
| Logarithmic $I_\nu$ and $P_\nu$ | $(@bind log_dust_positive PlutoUI.CheckBox(default = true)) |
| Symmetric-logarithmic $Q_\nu$ and $U_\nu$ | $(@bind log_dust_signed PlutoUI.CheckBox(default = false)) |
| Polarization pseudo-vectors over $I_\nu$ | $(@bind show_dust_vectors PlutoUI.CheckBox(default = true)) |
| Pseudo-vector stride | $(@bind dust_vector_stride PlutoUI.Slider(2:1:16; default = 5, show_value = true)) |
| $p$ versus $N_{\mathrm H}$ relation | $(@bind show_dust_p_column PlutoUI.CheckBox(default = true)) |
| Polarization-fraction PDF | $(@bind show_dust_p_pdf PlutoUI.CheckBox(default = true)) |

### Stokes spectrum at one sky pixel

The selected pixel is marked by a white cross on every dust map. The spectral panels evaluate the optically thin modified-blackbody model over the requested frequency interval. A vertical line and marker identify the **observing frequency** selected above. With the present single-temperature, single-$\beta_{\mathrm d}$ prescription, $I_\nu$, $Q_\nu$, and $U_\nu$ share the same frequency scaling while retaining their line-of-sight geometric amplitudes and signs.

| Pixel-spectrum setting | Control |
|:--|:--|
| First sky-axis pixel | $(@bind dust_sky_i PlutoUI.Slider(1:size(cube.rho, sky_dims[1]); default = cld(size(cube.rho, sky_dims[1]), 2), show_value = true)) |
| Second sky-axis pixel | $(@bind dust_sky_j PlutoUI.Slider(1:size(cube.rho, sky_dims[2]); default = cld(size(cube.rho, sky_dims[2]), 2), show_value = true)) |
| Minimum spectral frequency [$\mathrm{GHz}$] | $(@bind dust_spectrum_min_GHz PlutoUI.NumberField(1.0:1.0:2000.0; default = 30.0)) |
| Maximum spectral frequency [$\mathrm{GHz}$] | $(@bind dust_spectrum_max_GHz PlutoUI.NumberField(2.0:1.0:3000.0; default = 1200.0)) |
| Number of frequency samples | $(@bind dust_spectrum_samples PlutoUI.Select([64, 128, 256, 512]; default = 256)) |
| Logarithmic frequency axis | $(@bind log_dust_frequency_axis PlutoUI.CheckBox(default = true)) |
| $I_\nu$ spectrum | $(@bind show_dust_I_spectrum PlutoUI.CheckBox(default = true)) |
| $Q_\nu$ spectrum | $(@bind show_dust_Q_spectrum PlutoUI.CheckBox(default = true)) |
| $U_\nu$ spectrum | $(@bind show_dust_U_spectrum PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 130ccf03-d7cc-4b71-9210-dbc0e43dfa82
begin
    dust_nu_Hz = Float64(dust_frequency_GHz) * 1.0e9
    dust_T_K = max(Float64(dust_temperature_K), 2.8)
    dust_planck_MJysr = (2H_PLANCK_CGS * dust_nu_Hz^3 / C_LIGHT_CGS^2) /
        expm1(H_PLANCK_CGS * dust_nu_Hz / (K_B_CGS * dust_T_K)) * 1.0e17
    dust_sigma_cm2 = Float64(dust_sigma353_cm2) *
        (dust_nu_Hz / 353.0e9)^Float64(dust_beta)
    dust_nH = cube.rho ./ (max(Float64(dust_mu_H), eps(Float64)) * M_H_CGS)
    dust_Blos = Bcomponents[los_dim]
    dust_B1 = Bcomponents[sky_dims[1]]
    dust_B2 = Bcomponents[sky_dims[2]]
    dust_Bnorm2 = max.(mag.B2, eps(Float64))
    dust_cos2gamma = clamp.((dust_B1 .^ 2 .+ dust_B2 .^ 2) ./ dust_Bnorm2, 0.0, 1.0)
    dust_phi = atan.(dust_B2, dust_B1) .+ pi / 2
    dust_p0_value = Float64(dust_p0)
    dust_column_weight = dx_los_cm .* dust_nH
    dust_I_geometry = finite_sum_dims(dust_column_weight .*
        (1 .- dust_p0_value .* (dust_cos2gamma .- 2 / 3)),
        los_dim)
    dust_Q_geometry = finite_sum_dims(dust_column_weight .* dust_p0_value .* dust_cos2gamma .*
        cos.(2 .* dust_phi), los_dim)
    dust_U_geometry = finite_sum_dims(dust_column_weight .* dust_p0_value .* dust_cos2gamma .*
        sin.(2 .* dust_phi), los_dim)
    dust_I_geometry = apply_observational_beam_2d(dust_I_geometry, cube, sky_dims)
    dust_Q_geometry = apply_observational_beam_2d(dust_Q_geometry, cube, sky_dims)
    dust_U_geometry = apply_observational_beam_2d(dust_U_geometry, cube, sky_dims)
    dust_spectral_factor = dust_planck_MJysr * dust_sigma_cm2
    dust_I = dust_spectral_factor .* dust_I_geometry
    dust_Q = dust_spectral_factor .* dust_Q_geometry
    dust_U = dust_spectral_factor .* dust_U_geometry
    dust_P = sqrt.(dust_Q .^ 2 .+ dust_U .^ 2)
    dust_fraction = dust_P ./ max.(dust_I, eps(Float64))
    dust_angle_deg = rad2deg.(0.5 .* atan.(dust_U, dust_Q))

    dust_map_specs = NamedTuple[]
    show_dust_I && push!(dust_map_specs, (data = log_dust_positive ? safe_log10.(dust_I) : dust_I,
        label = log_dust_positive ? L"\log_{10}(I_\nu/[\mathrm{MJy\,sr}^{-1}])" : L"I_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        colormap = :magma, diverging = false, fixed_range = nothing, vectors = show_dust_vectors))
    show_dust_Q && push!(dust_map_specs, (data = log_dust_signed ? symlog10(dust_Q) : dust_Q,
        label = log_dust_signed ? L"\operatorname{symlog}_{10}(Q_\nu/[\mathrm{MJy\,sr}^{-1}])" : L"Q_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        colormap = :balance, diverging = true, fixed_range = nothing, vectors = false))
    show_dust_U && push!(dust_map_specs, (data = log_dust_signed ? symlog10(dust_U) : dust_U,
        label = log_dust_signed ? L"\operatorname{symlog}_{10}(U_\nu/[\mathrm{MJy\,sr}^{-1}])" : L"U_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        colormap = :balance, diverging = true, fixed_range = nothing, vectors = false))
    show_dust_P && push!(dust_map_specs, (data = log_dust_positive ? safe_log10.(dust_P) : dust_P,
        label = log_dust_positive ? L"\log_{10}(P_\nu/[\mathrm{MJy\,sr}^{-1}])" : L"P_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        colormap = :viridis, diverging = false, fixed_range = nothing, vectors = false))
    show_dust_fraction && push!(dust_map_specs, (data = dust_fraction,
        label = L"p=P_\nu/I_\nu", colormap = :plasma, diverging = false,
        fixed_range = (0.0, max(finite_quantile(dust_fraction, 0.995; default = 1.0e-6),
            1.0e-6)), vectors = false))
    show_dust_angle && push!(dust_map_specs, (data = dust_angle_deg,
        label = L"\psi\;[{}^\circ]", colormap = :hsv, diverging = false,
        fixed_range = (-90.0, 90.0), vectors = false))

    if isempty(dust_map_specs)
        fig_dust = Figure(size = (900, 180))
        Label(fig_dust[1, 1], L"\mathrm{Select\ at\ least\ one\ dust\ product.}", fontsize = 20)
    else
        dust_ncols = length(dust_map_specs) == 1 ? 1 : 2
        dust_nrows = cld(length(dust_map_specs), dust_ncols)
        fig_dust = Figure(size = (560dust_ncols, 420dust_nrows))
        for (index, spec) in enumerate(dust_map_specs)
            row, col = cld(index, dust_ncols), mod1(index, dust_ncols)
            panel = fig_dust[row, col] = GridLayout()
            ax = latex_axis(panel[1, 1], xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            dust_range = isnothing(spec.fixed_range) ?
                robust_colorrange(spec.data, color_percentile; diverging = spec.diverging) : spec.fixed_range
            hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap, colorrange = dust_range)
            scatter!(ax,
                [sky_coordinates[1][Int(dust_sky_i)]],
                [sky_coordinates[2][Int(dust_sky_j)]];
                marker = :cross, markersize = 18, strokewidth = 3,
                color = :white)
            latex_colorbar(panel[1, 2], hm; label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
            if spec.vectors
                stride = Int(dust_vector_stride)
                ix = collect(1:stride:size(dust_angle_deg, 1))
                iy = collect(1:stride:size(dust_angle_deg, 2))
                points = [Point2f(sky_coordinates[1][i], sky_coordinates[2][j]) for i in ix for j in iy]
                directions = [Vec2f(cosd(dust_angle_deg[i, j]), sind(dust_angle_deg[i, j])) for i in ix for j in iy]
                arrows2d!(ax, points, directions; normalize = true,
                    lengthscale = 0.65stride * step(sky_coordinates[1]), align = :center,
                    color = (:white, 0.88), shaftwidth = 1.4, tipwidth = 0, tiplength = 0)
            end
        end
    end
    display_dust_maps ? fig_dust : nothing
end

# ╔═╡ 5b3f6a91-246e-4cc3-8f68-164f7ff2f07c
begin
    dust_selected_frequency_GHz = Float64(dust_frequency_GHz)
    dust_frequency_low_GHz = max(min(Float64(dust_spectrum_min_GHz),
        Float64(dust_spectrum_max_GHz), dust_selected_frequency_GHz), eps(Float64))
    dust_frequency_high_GHz = max(max(Float64(dust_spectrum_min_GHz),
        Float64(dust_spectrum_max_GHz), dust_selected_frequency_GHz),
        dust_frequency_low_GHz * (1 + sqrt(eps(Float64))))
    dust_frequency_axis_GHz = log_dust_frequency_axis ?
        10.0 .^ range(log10(dust_frequency_low_GHz), log10(dust_frequency_high_GHz);
            length = Int(dust_spectrum_samples)) :
        collect(range(dust_frequency_low_GHz, dust_frequency_high_GHz;
            length = Int(dust_spectrum_samples)))
    dust_frequency_axis_Hz = dust_frequency_axis_GHz .* 1.0e9
    dust_planck_spectrum_MJysr = @. (2H_PLANCK_CGS * dust_frequency_axis_Hz^3 /
        C_LIGHT_CGS^2) / expm1(H_PLANCK_CGS * dust_frequency_axis_Hz /
        (K_B_CGS * dust_T_K)) * 1.0e17
    dust_sigma_spectrum_cm2 = @. Float64(dust_sigma353_cm2) *
        (dust_frequency_axis_Hz / 353.0e9)^Float64(dust_beta)
    dust_pixel_spectral_factor = dust_planck_spectrum_MJysr .* dust_sigma_spectrum_cm2
    dust_pixel_index = (Int(dust_sky_i), Int(dust_sky_j))
    dust_pixel_I_spectrum = dust_pixel_spectral_factor .* dust_I_geometry[dust_pixel_index...]
    dust_pixel_Q_spectrum = dust_pixel_spectral_factor .* dust_Q_geometry[dust_pixel_index...]
    dust_pixel_U_spectrum = dust_pixel_spectral_factor .* dust_U_geometry[dust_pixel_index...]

    dust_pixel_spectrum_specs = NamedTuple[]
    show_dust_I_spectrum && push!(dust_pixel_spectrum_specs, (
        values = dust_pixel_I_spectrum,
        selected_value = dust_I[dust_pixel_index...],
        ylabel = L"I_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        color = MHD_COLORS[1], signed = false))
    show_dust_Q_spectrum && push!(dust_pixel_spectrum_specs, (
        values = dust_pixel_Q_spectrum,
        selected_value = dust_Q[dust_pixel_index...],
        ylabel = L"Q_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        color = MHD_COLORS[2], signed = true))
    show_dust_U_spectrum && push!(dust_pixel_spectrum_specs, (
        values = dust_pixel_U_spectrum,
        selected_value = dust_U[dust_pixel_index...],
        ylabel = L"U_\nu\;[\mathrm{MJy\,sr}^{-1}]",
        color = MHD_COLORS[3], signed = true))

    if isempty(dust_pixel_spectrum_specs)
        fig_dust_pixel_spectrum = Figure(size = (900, 180))
        Label(fig_dust_pixel_spectrum[1, 1],
            L"\mathrm{Select\ at\ least\ one\ dust\ Stokes\ spectrum.}", fontsize = 20)
    else
        fig_dust_pixel_spectrum = Figure(
            size = (430length(dust_pixel_spectrum_specs), 390))
        for (index, spec) in enumerate(dust_pixel_spectrum_specs)
            axis = latex_axis(fig_dust_pixel_spectrum[1, index],
                xlabel = L"\nu\;[\mathrm{GHz}]", ylabel = spec.ylabel,
                xscale = log_dust_frequency_axis ? log10 : identity,
                xticks = log_dust_frequency_axis ? DECADE_TICKS : Makie.automatic,
                xminorticks = log_dust_frequency_axis ? IntervalsBetween(9) : IntervalsBetween(5),
                xminorticksvisible = true)
            spec.signed && hlines!(axis, [0.0];
                color = (:gray45, 0.55), linestyle = :dash, linewidth = 1.4)
            lines!(axis, dust_frequency_axis_GHz, spec.values;
                color = spec.color, linewidth = 2.7)
            vlines!(axis, [dust_selected_frequency_GHz];
                color = :black, linestyle = :dot, linewidth = 1.8)
            scatter!(axis, [dust_selected_frequency_GHz], [spec.selected_value];
                color = spec.color, marker = :star5, markersize = 14,
                strokecolor = :black, strokewidth = 0.8)
        end
    end
    display_dust_pixel_spectrum ? fig_dust_pixel_spectrum : nothing
end

# ╔═╡ 61c3b28c-9d62-4689-a85d-bc827e89641d
begin
    function dust_distribution_products(c)
        local_mag = magnetic_fields(c)
        local_Bcomponents = (c.bx, c.by, c.bz)
        local_B1 = local_Bcomponents[sky_dims[1]]
        local_B2 = local_Bcomponents[sky_dims[2]]
        local_cos2gamma = clamp.((local_B1 .^ 2 .+ local_B2 .^ 2) ./
            max.(local_mag.B2, eps(Float64)), 0.0, 1.0)
        local_phi = atan.(local_B2, local_B1) .+ pi / 2
        local_dx_cm = c.L[los_dim] / size(c.rho, los_dim) * PC_CM
        local_nH = c.rho ./ (max(Float64(dust_mu_H), eps(Float64)) * M_H_CGS)
        local_weight = local_dx_cm .* local_nH
        local_p0 = Float64(dust_p0)
        local_I = finite_sum_dims(local_weight .*
            (1 .- local_p0 .* (local_cos2gamma .- 2 / 3)), los_dim)
        local_Q = finite_sum_dims(local_weight .* local_p0 .* local_cos2gamma .*
            cos.(2 .* local_phi), los_dim)
        local_U = finite_sum_dims(local_weight .* local_p0 .* local_cos2gamma .*
            sin.(2 .* local_phi), los_dim)
        local_I = apply_observational_beam_2d(local_I, c, sky_dims)
        local_Q = apply_observational_beam_2d(local_Q, c, sky_dims)
        local_U = apply_observational_beam_2d(local_U, c, sky_dims)
        local_fraction = sqrt.(local_Q .^ 2 .+ local_U .^ 2) ./
            max.(local_I, eps(Float64))
        local_column = finite_sum_dims(c.rho, los_dim) .* local_dx_cm ./
            (max(Float64(dust_mu_H), eps(Float64)) * M_H_CGS)
        (column = local_column, fraction = local_fraction)
    end

    dust_distributions_by_run = Dict(label =>
        dust_distribution_products(comparison_cube(label))
        for label in comparison_run_labels)
    dust_stat_specs = Symbol[]
    show_dust_p_column && push!(dust_stat_specs, :column_relation)
    show_dust_p_pdf && push!(dust_stat_specs, :fraction_pdf)
    if isempty(dust_stat_specs)
        fig_dust_statistics = Figure(size = (900, 180))
        Label(fig_dust_statistics[1, 1], L"\mathrm{Select\ at\ least\ one\ dust\ statistic.}", fontsize = 20)
    else
        fig_dust_statistics = Figure(size = (540length(dust_stat_specs), 470))
        dust_distribution_vectors = Dict(label => begin
            products = dust_distributions_by_run[label]
            local_N = vec(Float64.(products.column))
            local_p = 100 .* vec(Float64.(products.fraction))
            valid = isfinite.(local_N) .& isfinite.(local_p) .&
                (local_N .> 0) .& (local_p .>= 0)
            (N = local_N[valid], p = local_p[valid])
        end for label in comparison_run_labels)
        all_dust_p = vcat([dust_distribution_vectors[label].p
            for label in comparison_run_labels]...)
        dust_pdf_upper = isempty(all_dust_p) ? 1.0e-6 :
            max(quantile(all_dust_p, 0.999), 1.0e-6)
        dust_pdf_edges = range(0, dust_pdf_upper; length = 45)
        for (index, statistic) in enumerate(dust_stat_specs)
            if statistic == :column_relation
                ax = latex_axis(fig_dust_statistics[1, index], xlabel = L"N_{\mathrm H}\;[\mathrm{cm}^{-2}]",
                    ylabel = L"100p_{\mathrm d}\;[\%]", xscale = log10,
                    xticks = DECADE_TICKS,
                    xminorticks = IntervalsBetween(9), xminorticksvisible = true)
                for label in comparison_run_labels
                    dust_N_valid = dust_distribution_vectors[label].N
                    dust_p_valid = dust_distribution_vectors[label].p
                    sample_step = max(1, cld(length(dust_N_valid), 6000))
                    sample = 1:sample_step:length(dust_N_valid)
                    scatter!(ax, dust_N_valid[sample], dust_p_valid[sample];
                        color = (run_colors[label], 0.16), markersize = 4)
                    if !isempty(dust_N_valid)
                        edges = 10 .^ range(extrema(log10.(dust_N_valid))...; length = 18)
                        centers, medians = Float64[], Float64[]
                        for bin in 1:length(edges)-1
                            members = (dust_N_valid .>= edges[bin]) .&
                                (dust_N_valid .< edges[bin + 1])
                            any(members) || continue
                            push!(centers, sqrt(edges[bin] * edges[bin + 1]))
                            push!(medians, median(dust_p_valid[members]))
                        end
                        lines!(ax, centers, medians; color = run_colors[label],
                            linewidth = 3, label = legend_run_label(label))
                        scatter!(ax, centers, medians;
                            color = run_colors[label], markersize = 7)
                    end
                end
            else
                ax = latex_axis(fig_dust_statistics[1, index], xlabel = L"100p_{\mathrm d}\;[\%]",
                    ylabel = L"\mathrm{PDF}(100p_{\mathrm d})")
                for label in comparison_run_labels
                    dust_p_valid = dust_distribution_vectors[label].p
                    isempty(dust_p_valid) && continue
                    histogram = fit(Histogram,
                        clamp.(dust_p_valid, 0, prevfloat(dust_pdf_upper)),
                        dust_pdf_edges)
                    centers = (dust_pdf_edges[1:end-1] .+ dust_pdf_edges[2:end]) ./ 2
                    probability = histogram.weights ./
                        max(sum(histogram.weights .* diff(dust_pdf_edges)), eps(Float64))
                    stairs!(ax, centers, probability; color = run_colors[label],
                        linewidth = 2.5, step = :center,
                        label = legend_run_label(label))
                end
            end
        end
        Legend(fig_dust_statistics[2, 1:length(dust_stat_specs)],
            [LineElement(color = run_colors[label], linewidth = 2.5)
                for label in comparison_run_labels],
            legend_run_label.(comparison_run_labels), L"\mathrm{Simulation}";
            orientation = :horizontal, tellheight = true, framevisible = false)
    end
    fig_dust_p_column = polarization_column_figure(
        column_density, dust_fraction,
        L"100p_{\mathrm d}\;[\%]", MHD_COLORS[2])
    display_dust_statistics ? fig_dust_statistics : nothing
end

# ╔═╡ a0010001-6f8c-4d0c-9a10-000000000001
begin
    dust_structure_specs = [
        (data = dust_I, label = L"I_\nu\;[\mathrm{MJy\,sr}^{-1}]", color = MHD_COLORS[1], period = nothing),
        (data = dust_Q, label = L"Q_\nu\;[\mathrm{MJy\,sr}^{-1}]", color = MHD_COLORS[2], period = nothing),
        (data = dust_U, label = L"U_\nu\;[\mathrm{MJy\,sr}^{-1}]", color = MHD_COLORS[3], period = nothing),
        (data = dust_P, label = L"P_\nu\;[\mathrm{MJy\,sr}^{-1}]", color = MHD_COLORS[4], period = nothing),
        (data = dust_fraction, label = L"p_{\mathrm d}", color = MHD_COLORS[5], period = nothing),
        (data = dust_angle_deg, label = L"\psi_{\mathrm d}\;[{}^\circ]", color = MHD_COLORS[6], period = 180.0),
    ]
    fig_dust_structure = display_observational_structure_functions ?
        observational_structure_figure(dust_structure_specs, cube, sky_dims,
            observational_structure_order, observational_structure_samples;
            heading = "Dust observable structure functions") : Figure(size = (900, 120))
    display_observational_structure_functions ? fig_dust_structure : nothing
end

# ╔═╡ 6f4e2d11-2a88-41f4-93dc-01b51d86fb4f
md"""
---

## 14. Dichroic starlight polarization

Cell-by-cell dichroic Mueller propagation toward background stars at the selected distance.

| Starlight figure | Display |
|:--|:--:|
| Final Stokes and polarization maps | $(@bind display_starlight_maps PlutoUI.CheckBox(default = true)) |
| Selected sight-line profiles | $(@bind display_starlight_profiles PlutoUI.CheckBox(default = true)) |
| Polarization fraction versus $N_{\rm H}$ | $(@bind display_starlight_p_column PlutoUI.CheckBox(default = true)) |

| Physical setting | Control |
|:--|:--|
| Star distance [$\mathrm{pc}$] | $(@bind starlight_star_distance_pc PlutoUI.NumberField(default = cube.L[los_dim])) |
| Intrinsic dichroic polarization $p_0$ | $(@bind starlight_p0 PlutoUI.NumberField(0.0:0.005:0.95; default = 0.20)) |
| $N_{\rm H}/A_V$ [$\mathrm{cm}^{-2}$] | $(@bind starlight_nh_per_av PlutoUI.NumberField(default = 1.8e21)) |
| Gas mass per H nucleon [$m_{\rm H}$] | $(@bind starlight_mu_H PlutoUI.NumberField(1.0:0.01:2.0; default = 1.4)) |
| Incident $I_0$ | $(@bind starlight_initial_I PlutoUI.NumberField(default = 1.0)) |
| Incident $Q_0$ | $(@bind starlight_initial_Q PlutoUI.NumberField(default = 0.0)) |
| Incident $U_0$ | $(@bind starlight_initial_U PlutoUI.NumberField(default = 0.0)) |
| Incident $V_0$ | $(@bind starlight_initial_V PlutoUI.NumberField(default = 0.0)) |
| First sky-axis pixel | $(@bind starlight_sky_i PlutoUI.Slider(1:size(cube.rho, sky_dims[1]); default = cld(size(cube.rho, sky_dims[1]), 2), show_value = true)) |
| Second sky-axis pixel | $(@bind starlight_sky_j PlutoUI.Slider(1:size(cube.rho, sky_dims[2]); default = cld(size(cube.rho, sky_dims[2]), 2), show_value = true)) |

| Final map | Display |
|:--|:--:|
| Transmitted intensity $I/I_0$ | $(@bind show_starlight_I_map PlutoUI.CheckBox(default = true)) |
| Stokes $Q/I_0$ | $(@bind show_starlight_Q_map PlutoUI.CheckBox(default = false)) |
| Stokes $U/I_0$ | $(@bind show_starlight_U_map PlutoUI.CheckBox(default = false)) |
| Polarization fraction $p_\star$ | $(@bind show_starlight_p_map PlutoUI.CheckBox(default = true)) |
| Polarization angle $\psi_\star$ | $(@bind show_starlight_angle_map PlutoUI.CheckBox(default = true)) |
| Effective optical depth $\tau_V$ | $(@bind show_starlight_tau_map PlutoUI.CheckBox(default = true)) |
| Density-weighted $B_{\rm LOS}$ | $(@bind show_starlight_blos_map PlutoUI.CheckBox(default = false)) |

| Selected-pixel profile | Display |
|:--|:--:|
| $I/I_0$ | $(@bind show_starlight_I_profile PlutoUI.CheckBox(default = true)) |
| $Q/I_0$ | $(@bind show_starlight_Q_profile PlutoUI.CheckBox(default = true)) |
| $U/I_0$ | $(@bind show_starlight_U_profile PlutoUI.CheckBox(default = true)) |
| $p_\star$ | $(@bind show_starlight_p_profile PlutoUI.CheckBox(default = true)) |
| $\psi_\star$ | $(@bind show_starlight_angle_profile PlutoUI.CheckBox(default = true)) |
| Cumulative $\tau_V$ | $(@bind show_starlight_tau_profile PlutoUI.CheckBox(default = true)) |
| Local $B_{\rm LOS}$ | $(@bind show_starlight_blos_profile PlutoUI.CheckBox(default = false)) |
| Local magnetic inclination $\gamma_B$ | $(@bind show_starlight_gamma_profile PlutoUI.CheckBox(default = false)) |
"""

# ╔═╡ a2d5319b-06f4-4efc-b3d7-3a9719292305
begin
    function starlight_mueller_step(I, Q, U, V, tau, psi, p0)
        delta_tau_factor = log((1 + p0) / (1 - p0))
        delta_tau = tau .* delta_tau_factor
        qtrans = exp.(-(tau .+ delta_tau))
        rtrans = exp.(-(tau .- delta_tau))
        sin2psi = sin.(2 .* psi)
        cos2psi = cos.(2 .* psi)
        qplus = qtrans .+ rtrans
        qminus = qtrans .- rtrans
        root = sqrt.(max.(qtrans .* rtrans, 0.0))
        cross = cos2psi .* sin2psi
        Inew = 0.5 .* (qplus .* I .+ qminus .* cos2psi .* Q .+
            qminus .* sin2psi .* U)
        Qnew = 0.5 .* (qminus .* cos2psi .* I .+
            (qplus .* cos2psi .^ 2 .+ 2 .* root .* sin2psi .^ 2) .* Q .+
            (qplus .- 2 .* root) .* cross .* U)
        Unew = 0.5 .* (qminus .* sin2psi .* I .+
            (qplus .- 2 .* root) .* cross .* Q .+
            (qplus .* sin2psi .^ 2 .+ 2 .* root .* cos2psi .^ 2) .* U)
        Vnew = root .* V
        Inew, Qnew, Unew, Vnew
    end

    starlight_nlos = size(cube.rho, los_dim)
    starlight_dx_pc = cube.L[los_dim] / starlight_nlos
    starlight_requested_distance_pc = clamp(Float64(starlight_star_distance_pc),
        starlight_dx_pc, cube.L[los_dim])
    starlight_cell_count = clamp(ceil(Int,
        starlight_requested_distance_pc / starlight_dx_pc), 1, starlight_nlos)
    starlight_distance_pc = collect(1:starlight_cell_count) .* starlight_dx_pc
    starlight_actual_distance_pc = last(starlight_distance_pc)
    starlight_shape = size(selectdim(cube.rho, los_dim, 1))
    starlight_I0 = Float64(starlight_initial_I)
    starlight_I0 > 0 || error("Incident starlight intensity I0 must be positive.")
    starlight_Q0 = Float64(starlight_initial_Q)
    starlight_U0 = Float64(starlight_initial_U)
    starlight_V0 = Float64(starlight_initial_V)
    starlight_p0_value = clamp(Float64(starlight_p0), 0.0, 1 - eps(Float64))
    starlight_nh_av = max(Float64(starlight_nh_per_av), eps(Float64))
    starlight_I_map = fill(starlight_I0, starlight_shape)
    starlight_Q_map = fill(starlight_Q0, starlight_shape)
    starlight_U_map = fill(starlight_U0, starlight_shape)
    starlight_V_map = fill(starlight_V0, starlight_shape)
    starlight_tau_map = zeros(Float64, starlight_shape)
    starlight_blos_numerator = zeros(Float64, starlight_shape)
    starlight_column = zeros(Float64, starlight_shape)
    starlight_pixel_index = (Int(starlight_sky_i), Int(starlight_sky_j))
    starlight_I_profile = Float64[]
    starlight_Q_profile = Float64[]
    starlight_U_profile = Float64[]
    starlight_V_profile = Float64[]
    starlight_tau_profile = Float64[]
    starlight_blos_profile = Float64[]
    starlight_gamma_profile_deg = Float64[]
    starlight_local_psi_deg = Float64[]
    starlight_east_component = Bcomponents[sky_dims[1]]
    starlight_north_component = Bcomponents[sky_dims[2]]

    for cell_index in 1:starlight_cell_count
        density_slice = selectdim(cube.rho, los_dim, cell_index) ./
            (max(Float64(starlight_mu_H), eps(Float64)) * M_H_CGS)
        blos_slice = selectdim(Bcomponents[los_dim], los_dim, cell_index)
        beast_slice = selectdim(starlight_east_component, los_dim, cell_index)
        bnorth_slice = selectdim(starlight_north_component, los_dim, cell_index)
        bnorm_slice = sqrt.(blos_slice .^ 2 .+ beast_slice .^ 2 .+ bnorth_slice .^ 2)
        psi_slice = atan.(beast_slice, bnorth_slice)
        tau_slice = density_slice .* (starlight_dx_pc * PC_CM) ./ starlight_nh_av
        next_I_map, next_Q_map, next_U_map, next_V_map =
            starlight_mueller_step(starlight_I_map, starlight_Q_map,
                starlight_U_map, starlight_V_map, tau_slice, psi_slice,
                starlight_p0_value)
        starlight_I_map .= next_I_map
        starlight_Q_map .= next_Q_map
        starlight_U_map .= next_U_map
        starlight_V_map .= next_V_map
        starlight_tau_map .+= tau_slice
        starlight_blos_numerator .+= density_slice .* blos_slice
        starlight_column .+= density_slice
        push!(starlight_I_profile, starlight_I_map[starlight_pixel_index...])
        push!(starlight_Q_profile, starlight_Q_map[starlight_pixel_index...])
        push!(starlight_U_profile, starlight_U_map[starlight_pixel_index...])
        push!(starlight_V_profile, starlight_V_map[starlight_pixel_index...])
        push!(starlight_tau_profile, starlight_tau_map[starlight_pixel_index...])
        push!(starlight_blos_profile,
            GAUSS_TO_MICROGAUSS * blos_slice[starlight_pixel_index...])
        local_B = bnorm_slice[starlight_pixel_index...]
        local_blos = blos_slice[starlight_pixel_index...]
        push!(starlight_gamma_profile_deg, local_B > 0 ?
            rad2deg(acos(clamp(local_blos / local_B, -1.0, 1.0))) : NaN)
        push!(starlight_local_psi_deg,
            rad2deg(psi_slice[starlight_pixel_index...]))
    end

    starlight_I_map = apply_observational_beam_2d(starlight_I_map, cube, sky_dims)
    starlight_Q_map = apply_observational_beam_2d(starlight_Q_map, cube, sky_dims)
    starlight_U_map = apply_observational_beam_2d(starlight_U_map, cube, sky_dims)
    starlight_V_map = apply_observational_beam_2d(starlight_V_map, cube, sky_dims)
    starlight_I_normalized = starlight_I_map ./ starlight_I0
    starlight_Q_normalized = starlight_Q_map ./ starlight_I0
    starlight_U_normalized = starlight_U_map ./ starlight_I0
    starlight_p_map = clamp.(sqrt.(starlight_Q_map .^ 2 .+
        starlight_U_map .^ 2) ./ max.(abs.(starlight_I_map), eps(Float64)), 0.0, 1.0)
    starlight_angle_deg = rad2deg.(0.5 .* atan.(starlight_U_map, starlight_Q_map))
    starlight_NH_map = starlight_column .* (starlight_dx_pc * PC_CM)
    starlight_blos_map_uG = GAUSS_TO_MICROGAUSS .* starlight_blos_numerator ./
        max.(starlight_column, eps(Float64))
    starlight_I_profile_normalized = starlight_I_profile ./ starlight_I0
    starlight_Q_profile_normalized = starlight_Q_profile ./ starlight_I0
    starlight_U_profile_normalized = starlight_U_profile ./ starlight_I0
    starlight_p_profile = clamp.(sqrt.(starlight_Q_profile .^ 2 .+
        starlight_U_profile .^ 2) ./ max.(abs.(starlight_I_profile), eps(Float64)), 0.0, 1.0)
    starlight_angle_profile_deg = rad2deg.(0.5 .* atan.(
        starlight_U_profile, starlight_Q_profile))
end

# ╔═╡ d1e7626c-81d7-48cd-90be-695a13aa9997
begin
    starlight_final_I = starlight_I_normalized[starlight_pixel_index...]
    starlight_final_Q = starlight_Q_normalized[starlight_pixel_index...]
    starlight_final_U = starlight_U_normalized[starlight_pixel_index...]
    starlight_final_p = starlight_p_map[starlight_pixel_index...]
    starlight_final_angle = starlight_angle_deg[starlight_pixel_index...]
    starlight_final_tau = starlight_tau_map[starlight_pixel_index...]
    Markdown.parse("""
    ### Active starlight sight line

    | Quantity | Value |
    |:--|:--|
    | Sky pixel | **($(starlight_sky_i), $(starlight_sky_j))** |
    | Integrated stellar distance | **$(@sprintf("%.6g", starlight_actual_distance_pc))** ``\\mathrm{pc}`` |
    | Final ``I/I_0`` | **$(@sprintf("%.6g", starlight_final_I))** |
    | Final ``Q/I_0`` | **$(@sprintf("%.6g", starlight_final_Q))** |
    | Final ``U/I_0`` | **$(@sprintf("%.6g", starlight_final_U))** |
    | Polarization fraction ``p_\\star`` | **$(@sprintf("%.6g", starlight_final_p))** |
    | Polarization angle ``\\psi_\\star`` | **$(@sprintf("%.6g", starlight_final_angle))** ``{}^\\circ`` |
    | Effective optical depth ``\\tau_V`` | **$(@sprintf("%.6g", starlight_final_tau))** |
    """)
end

# ╔═╡ bcc9a131-9df2-455a-bdac-586323fd58d0
begin
    starlight_map_specs = NamedTuple[]
    show_starlight_I_map && push!(starlight_map_specs, (
        data = starlight_I_normalized, label = L"I/I_0", colormap = :magma,
        diverging = false, fixed_range = nothing))
    show_starlight_Q_map && push!(starlight_map_specs, (
        data = starlight_Q_normalized, label = L"Q/I_0", colormap = :balance,
        diverging = true, fixed_range = nothing))
    show_starlight_U_map && push!(starlight_map_specs, (
        data = starlight_U_normalized, label = L"U/I_0", colormap = :balance,
        diverging = true, fixed_range = nothing))
    show_starlight_p_map && push!(starlight_map_specs, (
        data = starlight_p_map, label = L"p_\star", colormap = :plasma,
        diverging = false, fixed_range = (0.0,
            max(finite_quantile(starlight_p_map, 0.995; default = 1.0e-6), 1.0e-6))))
    show_starlight_angle_map && push!(starlight_map_specs, (
        data = starlight_angle_deg, label = L"\psi_\star\;[{}^\circ]", colormap = :hsv,
        diverging = false, fixed_range = (-90.0, 90.0)))
    show_starlight_tau_map && push!(starlight_map_specs, (
        data = starlight_tau_map, label = L"\tau_V", colormap = :viridis,
        diverging = false, fixed_range = nothing))
    show_starlight_blos_map && push!(starlight_map_specs, (
        data = starlight_blos_map_uG,
        label = L"\langle B_{\mathrm{LOS}}\rangle_n\;[\mu\mathrm{G}]",
        colormap = :balance, diverging = true, fixed_range = nothing))

    if isempty(starlight_map_specs)
        fig_starlight_maps = Figure(size = (900, 180))
        Label(fig_starlight_maps[1, 1],
            L"\mathrm{Select\ at\ least\ one\ starlight\ map.}", fontsize = 20)
    else
        starlight_map_ncols = length(starlight_map_specs) == 1 ? 1 : 2
        starlight_map_nrows = cld(length(starlight_map_specs), starlight_map_ncols)
        fig_starlight_maps = Figure(
            size = (560starlight_map_ncols, 420starlight_map_nrows))
        for (index, spec) in enumerate(starlight_map_specs)
            row, col = cld(index, starlight_map_ncols), mod1(index, starlight_map_ncols)
            panel = fig_starlight_maps[row, col] = GridLayout()
            axis = latex_axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            colorrange = isnothing(spec.fixed_range) ?
                robust_colorrange(spec.data, color_percentile;
                    diverging = spec.diverging) : spec.fixed_range
            heat = heatmap!(axis, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap, colorrange)
            scatter!(axis,
                [sky_coordinates[1][Int(starlight_sky_i)]],
                [sky_coordinates[2][Int(starlight_sky_j)]];
                marker = :star5, markersize = 17, color = :white,
                strokecolor = :black, strokewidth = 1.2)
            latex_colorbar(panel[1, 2], heat; label = as_latex(spec.label),
                tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_starlight_maps ? fig_starlight_maps : nothing
end

# ╔═╡ 0c8f7453-06e8-47ea-ae0f-4f09a89b16ae
begin
    starlight_profile_specs = NamedTuple[]
    show_starlight_I_profile && push!(starlight_profile_specs, (
        values = starlight_I_profile_normalized, ylabel = L"I/I_0",
        color = MHD_COLORS[1], signed = false))
    show_starlight_Q_profile && push!(starlight_profile_specs, (
        values = starlight_Q_profile_normalized, ylabel = L"Q/I_0",
        color = MHD_COLORS[2], signed = true))
    show_starlight_U_profile && push!(starlight_profile_specs, (
        values = starlight_U_profile_normalized, ylabel = L"U/I_0",
        color = MHD_COLORS[3], signed = true))
    show_starlight_p_profile && push!(starlight_profile_specs, (
        values = starlight_p_profile, ylabel = L"p_\star",
        color = MHD_COLORS[4], signed = false))
    show_starlight_angle_profile && push!(starlight_profile_specs, (
        values = starlight_angle_profile_deg, ylabel = L"\psi_\star\;[{}^\circ]",
        color = MHD_COLORS[5], signed = true))
    show_starlight_tau_profile && push!(starlight_profile_specs, (
        values = starlight_tau_profile, ylabel = L"\tau_V(<d)",
        color = MHD_COLORS[6], signed = false))
    show_starlight_blos_profile && push!(starlight_profile_specs, (
        values = starlight_blos_profile,
        ylabel = L"B_{\mathrm{LOS}}\;[\mu\mathrm{G}]",
        color = MHD_COLORS[1], signed = true))
    show_starlight_gamma_profile && push!(starlight_profile_specs, (
        values = starlight_gamma_profile_deg,
        ylabel = L"\gamma_B\;[{}^\circ]",
        color = MHD_COLORS[4], signed = false))

    if isempty(starlight_profile_specs)
        fig_starlight_profiles = Figure(size = (900, 180))
        Label(fig_starlight_profiles[1, 1],
            L"\mathrm{Select\ at\ least\ one\ starlight\ profile.}", fontsize = 20)
    else
        starlight_profile_ncols = length(starlight_profile_specs) == 1 ? 1 : 2
        starlight_profile_nrows = cld(length(starlight_profile_specs),
            starlight_profile_ncols)
        fig_starlight_profiles = Figure(
            size = (540starlight_profile_ncols, 350starlight_profile_nrows))
        for (index, spec) in enumerate(starlight_profile_specs)
            row, col = cld(index, starlight_profile_ncols),
                mod1(index, starlight_profile_ncols)
            axis = latex_axis(fig_starlight_profiles[row, col],
                xlabel = L"d\;[\mathrm{pc}]", ylabel = spec.ylabel)
            spec.signed && hlines!(axis, [0.0]; color = (:gray45, 0.55),
                linestyle = :dash, linewidth = 1.3)
            valid = isfinite.(starlight_distance_pc) .& isfinite.(spec.values)
            lines!(axis, starlight_distance_pc[valid], spec.values[valid];
                color = spec.color, linewidth = 2.6)
            any(valid) && scatter!(axis,
                [starlight_distance_pc[findlast(valid)]],
                [spec.values[findlast(valid)]];
                color = spec.color, marker = :star5, markersize = 12,
                strokecolor = :black, strokewidth = 0.7)
        end
    end
    display_starlight_profiles ? fig_starlight_profiles : nothing
end

# ╔═╡ 7e8d4ac1-7b6e-44a8-981b-9c3a19c8de10
begin
    fig_starlight_p_column = polarization_column_figure(
        starlight_NH_map, starlight_p_map,
        L"100p_\star\;[\%]", MHD_COLORS[4])
    display_starlight_p_column ? fig_starlight_p_column : nothing
end

# ╔═╡ a0020002-6f8c-4d0c-9a10-000000000002
begin
    starlight_structure_specs = [
        (data = starlight_I_normalized, label = L"I/I_0", color = MHD_COLORS[1], period = nothing),
        (data = starlight_Q_normalized, label = L"Q/I_0", color = MHD_COLORS[2], period = nothing),
        (data = starlight_U_normalized, label = L"U/I_0", color = MHD_COLORS[3], period = nothing),
        (data = starlight_V_map ./ starlight_I0, label = L"V/I_0", color = MHD_COLORS[4], period = nothing),
        (data = starlight_p_map, label = L"p_\star", color = MHD_COLORS[5], period = nothing),
        (data = starlight_angle_deg, label = L"\psi_\star\;[{}^\circ]", color = MHD_COLORS[6], period = 180.0),
        (data = starlight_tau_map, label = L"\tau_V", color = MHD_COLORS[1], period = nothing),
        (data = starlight_blos_map_uG, label = L"\langle B_{\mathrm{LOS}}\rangle_n\;[\mu\mathrm{G}]", color = MHD_COLORS[2], period = nothing),
    ]
    fig_starlight_structure = display_observational_structure_functions ?
        observational_structure_figure(starlight_structure_specs, cube, sky_dims,
            observational_structure_order, observational_structure_samples;
            heading = "Starlight observable structure functions") : Figure(size = (900, 120))
    display_observational_structure_functions ? fig_starlight_structure : nothing
end

# ╔═╡ 67f95c39-1888-4d23-a2c2-2ee3a6cd7f0f
md"""
---

## 15. H I Zeeman splitting

Weak-splitting $\mathrm{H\,I}$ Zeeman synthesis. Select a sky pixel to compare the true weighted field with the Stokes-$V$ fit.

| Zeeman figure | Display |
|:--|:--:|
| Zeeman maps | $(@bind display_zeeman_maps PlutoUI.CheckBox(default = true)) |
| Stokes spectra | $(@bind display_zeeman_spectra PlutoUI.CheckBox(default = true)) |
| Circular-polarization fraction versus $N_{\rm HI}$ | $(@bind display_zeeman_p_column PlutoUI.CheckBox(default = true)) |

| Zeeman setting | Control |
|:--|:--|
| Neutral $\mathrm{H\,I}$ fraction | $(@bind zeeman_neutral_fraction PlutoUI.NumberField(0.0:0.01:1.0; default = 1.0)) |
| Rest frequency [$\mathrm{MHz}$] | $(@bind zeeman_frequency_MHz PlutoUI.NumberField(default = 1420.40575177)) |
| Splitting coefficient [$\mathrm{Hz}\,\mu\mathrm{G}^{-1}$] | $(@bind zeeman_coefficient_Hz_uG PlutoUI.NumberField(0.0:0.01:10.0; default = 2.80)) |
| Non-thermal line width [$\mathrm{km\,s}^{-1}$] | $(@bind zeeman_microturbulence_kms PlutoUI.NumberField(0.05:0.05:20.0; default = 0.8)) |
| Velocity padding [$\mathrm{km\,s}^{-1}$] | $(@bind zeeman_velocity_padding_kms PlutoUI.NumberField(1.0:0.5:50.0; default = 5.0)) |
| Number of channels | $(@bind zeeman_channel_count PlutoUI.Select([101, 201, 301, 401]; default = 201)) |
| First sky-axis pixel | $(@bind zeeman_sky_i PlutoUI.Slider(1:size(cube.rho, sky_dims[1]); default = cld(size(cube.rho, sky_dims[1]), 2), show_value = true)) |
| Second sky-axis pixel | $(@bind zeeman_sky_j PlutoUI.Slider(1:size(cube.rho, sky_dims[2]); default = cld(size(cube.rho, sky_dims[2]), 2), show_value = true)) |
| $\mathrm{H\,I}$-weighted $B_{\mathrm{LOS}}$ map | $(@bind show_zeeman_Bmap PlutoUI.CheckBox(default = true)) |
| Frequency-splitting map | $(@bind show_zeeman_split_map PlutoUI.CheckBox(default = true)) |
| Stokes-$I$ spectrum | $(@bind show_zeeman_I_spectrum PlutoUI.CheckBox(default = true)) |
| Stokes-$V$ spectrum | $(@bind show_zeeman_V_spectrum PlutoUI.CheckBox(default = true)) |
| Derivative-fit model | $(@bind show_zeeman_fit PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 76b9cf43-a13e-44c8-97d6-4760cc3aa486
begin
    function extract_sightline(A, i, j, line_dim, plane_dims)
        indices = Any[Colon(), Colon(), Colon()]
        indices[plane_dims[1]] = i
        indices[plane_dims[2]] = j
        vec(view(A, indices...))
    end

    zeeman_nHI = Float64(zeeman_neutral_fraction) .* cube.rho ./ (1.4M_H_CGS)
    zeeman_Blos_uG_cube = GAUSS_TO_MICROGAUSS .* Bcomponents[los_dim]
    zeeman_valid_cube = isfinite.(zeeman_nHI) .& isfinite.(zeeman_Blos_uG_cube) .&
        (zeeman_nHI .> 0)
    zeeman_weight_sum = finite_sum_dims(ifelse.(zeeman_valid_cube, zeeman_nHI, NaN), los_dim)
    zeeman_B_numerator = finite_sum_dims(ifelse.(zeeman_valid_cube,
        zeeman_nHI .* zeeman_Blos_uG_cube, NaN), los_dim)
    zeeman_weight_sum = apply_observational_beam_2d(zeeman_weight_sum, cube, sky_dims)
    zeeman_B_numerator = apply_observational_beam_2d(
        zeeman_B_numerator, cube, sky_dims)
    zeeman_Bmap_uG = zeeman_B_numerator ./ max.(zeeman_weight_sum, eps(Float64))
    zeeman_z = Float64(zeeman_coefficient_Hz_uG)
    zeeman_split_map_Hz = zeeman_z .* zeeman_Bmap_uG

    zeeman_i, zeeman_j = Int(zeeman_sky_i), Int(zeeman_sky_j)
    zeeman_v_components = (cube.vx, cube.vy, cube.vz)
    zeeman_vlos = extract_sightline(zeeman_v_components[los_dim], zeeman_i, zeeman_j, los_dim, sky_dims)
    zeeman_Tline = extract_sightline(T, zeeman_i, zeeman_j, los_dim, sky_dims)
    zeeman_nline = extract_sightline(zeeman_nHI, zeeman_i, zeeman_j, los_dim, sky_dims)
    zeeman_Bline_uG = extract_sightline(zeeman_Blos_uG_cube, zeeman_i, zeeman_j, los_dim, sky_dims)
    zeeman_sigma_thermal = sqrt.(K_B_CGS .* max.(zeeman_Tline, 0.0) ./ M_H_CGS) ./ KM_CM
    zeeman_sigma_kms = sqrt.(zeeman_sigma_thermal .^ 2 .+ Float64(zeeman_microturbulence_kms)^2)
    zeeman_Ncell = zeeman_nline .* dx_los_cm
    zeeman_line_valid = isfinite.(zeeman_vlos) .& isfinite.(zeeman_sigma_kms) .&
        isfinite.(zeeman_Ncell) .& isfinite.(zeeman_Bline_uG) .&
        (zeeman_sigma_kms .> 0) .& (zeeman_Ncell .>= 0)
    zeeman_vlos = zeeman_vlos[zeeman_line_valid]
    zeeman_sigma_kms = zeeman_sigma_kms[zeeman_line_valid]
    zeeman_Ncell = zeeman_Ncell[zeeman_line_valid]
    zeeman_Bline_uG = zeeman_Bline_uG[zeeman_line_valid]
    if isempty(zeeman_vlos)
        zeeman_vlos = [0.0]
        zeeman_sigma_kms = [max(Float64(zeeman_microturbulence_kms), eps(Float64))]
        zeeman_Ncell = [0.0]
        zeeman_Bline_uG = [0.0]
    end
    zeeman_padding = Float64(zeeman_velocity_padding_kms)
    zeeman_vmin = minimum(zeeman_vlos .- 4 .* zeeman_sigma_kms) - zeeman_padding
    zeeman_vmax = maximum(zeeman_vlos .+ 4 .* zeeman_sigma_kms) + zeeman_padding
    zeeman_velocity_axis = collect(range(zeeman_vmin, zeeman_vmax; length = Int(zeeman_channel_count)))
    zeeman_amplitudes = zeeman_Ncell ./
        (1.823e18 .* sqrt(2pi) .* max.(zeeman_sigma_kms, eps(Float64)))
    zeeman_profiles = [zeeman_amplitudes[cell] *
        exp(-0.5 * ((velocity - zeeman_vlos[cell]) / zeeman_sigma_kms[cell])^2)
        for velocity in zeeman_velocity_axis, cell in eachindex(zeeman_vlos)]
    zeeman_I_K = vec(sum(zeeman_profiles; dims = 2))
    zeeman_nu0_Hz = Float64(zeeman_frequency_MHz) * 1.0e6
    zeeman_dIdnu_profiles = [zeeman_profiles[channel, cell] *
        (zeeman_velocity_axis[channel] - zeeman_vlos[cell]) / zeeman_sigma_kms[cell]^2 *
        C_LIGHT_KMS / zeeman_nu0_Hz
        for channel in eachindex(zeeman_velocity_axis), cell in eachindex(zeeman_vlos)]
    zeeman_dIdnu = vec(sum(zeeman_dIdnu_profiles; dims = 2))
    zeeman_V_K = 0.5zeeman_z .* vec(zeeman_dIdnu_profiles * zeeman_Bline_uG)
    zeeman_fit_denominator = zeeman_z * sum(abs2, zeeman_dIdnu)
    zeeman_Bfit_uG = zeeman_fit_denominator > 0 ?
        2sum(zeeman_dIdnu .* zeeman_V_K) / zeeman_fit_denominator : NaN
    zeeman_Vfit_K = 0.5zeeman_z * zeeman_Bfit_uG .* zeeman_dIdnu
    zeeman_Btrue_uG = sum(zeeman_Ncell .* zeeman_Bline_uG) / max(sum(zeeman_Ncell), eps(Float64))

    zeeman_map_specs = NamedTuple[]
    show_zeeman_Bmap && push!(zeeman_map_specs, (data = zeeman_Bmap_uG,
        label = L"\langle B_{\mathrm{LOS}}\rangle_{\mathrm{HI}}\;[\mu\mathrm{G}]", colormap = :balance))
    show_zeeman_split_map && push!(zeeman_map_specs, (data = zeeman_split_map_Hz,
        label = L"\Delta\nu_Z\;[\mathrm{Hz}]", colormap = :balance))
    if isempty(zeeman_map_specs)
        fig_zeeman_maps = Figure(size = (900, 180))
        Label(fig_zeeman_maps[1, 1], L"\mathrm{Select\ at\ least\ one\ Zeeman\ map.}", fontsize = 20)
    else
        fig_zeeman_maps = Figure(size = (540length(zeeman_map_specs), 410))
        for (index, spec) in enumerate(zeeman_map_specs)
            panel = fig_zeeman_maps[1, index] = GridLayout()
            ax = latex_axis(panel[1, 1], xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap,
                colorrange = robust_colorrange(spec.data, color_percentile; diverging = true))
            scatter!(ax, [sky_coordinates[1][zeeman_i]], [sky_coordinates[2][zeeman_j]];
                marker = :cross, markersize = 20, strokewidth = 3, color = :white)
            latex_colorbar(panel[1, 2], hm; label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_zeeman_maps ? fig_zeeman_maps : nothing
end

# ╔═╡ 82e22c29-1cf5-47dc-a54c-68577f8069bc
begin
    zeeman_spectrum_specs = Symbol[]
    show_zeeman_I_spectrum && push!(zeeman_spectrum_specs, :I)
    show_zeeman_V_spectrum && push!(zeeman_spectrum_specs, :V)
    if isempty(zeeman_spectrum_specs)
        fig_zeeman_spectra = Figure(size = (900, 180))
        Label(fig_zeeman_spectra[1, 1], L"\mathrm{Select\ at\ least\ one\ Zeeman\ spectrum.}", fontsize = 20)
    else
        fig_zeeman_spectra = Figure(size = (560length(zeeman_spectrum_specs), 390))
        for (index, spectrum) in enumerate(zeeman_spectrum_specs)
            if spectrum == :I
                ax = latex_axis(fig_zeeman_spectra[1, index], xlabel = L"v_{\mathrm{LOS}}\;[\mathrm{km\,s}^{-1}]",
                    ylabel = L"I_\nu\;[\mathrm{K}]")
                lines!(ax, zeeman_velocity_axis, zeeman_I_K; color = MHD_COLORS[1], linewidth = 2.5)
            else
                ax = latex_axis(fig_zeeman_spectra[1, index], xlabel = L"v_{\mathrm{LOS}}\;[\mathrm{km\,s}^{-1}]",
                    ylabel = L"V_\nu\;[\mathrm{mK}]")
                lines!(ax, zeeman_velocity_axis, 1.0e3 .* zeeman_V_K; color = MHD_COLORS[2], linewidth = 2.5)
                show_zeeman_fit && lines!(ax, zeeman_velocity_axis, 1.0e3 .* zeeman_Vfit_K;
                    color = :black, linewidth = 2, linestyle = :dash)
                hlines!(ax, [0.0]; color = (:gray, 0.5), linestyle = :dot)
            end
        end
    end
    display_zeeman_spectra ? fig_zeeman_spectra : nothing
end

# ╔═╡ 8f9e5bd2-8c7f-45b9-a92c-ad4b20d9ef21
begin
    function zeeman_fraction_maps(nHI, vlos_cube, temperature_cube, Blos_uG_cube,
            line_dim, plane_dims, dx_cm, microturbulence_kms,
            rest_frequency_Hz, splitting_coefficient, channel_count)
        map_shape = size(selectdim(nHI, line_dim, 1))
        NHI_map = zeros(Float64, map_shape)
        pV_map = zeros(Float64, map_shape)
        Ipeak_map = zeros(Float64, map_shape)
        for pixel in CartesianIndices(map_shape)
            i, j = Tuple(pixel)
            nline = extract_sightline(nHI, i, j, line_dim, plane_dims)
            vline = extract_sightline(vlos_cube, i, j, line_dim, plane_dims)
            Tline = extract_sightline(temperature_cube, i, j, line_dim, plane_dims)
            Bline = extract_sightline(Blos_uG_cube, i, j, line_dim, plane_dims)
            sigma = sqrt.(K_B_CGS .* max.(Tline, 0.0) ./ M_H_CGS) ./ KM_CM
            sigma = sqrt.(sigma .^ 2 .+ microturbulence_kms^2)
            Ncell = nline .* dx_cm
            NHI_map[pixel] = sum(Ncell)
            amplitudes = Ncell ./
                (1.823e18 .* sqrt(2pi) .* max.(sigma, eps(Float64)))
            velocity = range(minimum(vline .- 5 .* sigma),
                maximum(vline .+ 5 .* sigma); length = channel_count)
            Iprofile = zeros(Float64, channel_count)
            Vprofile = zeros(Float64, channel_count)
            for cell in eachindex(vline)
                profile = @. amplitudes[cell] *
                    exp(-0.5 * ((velocity - vline[cell]) / sigma[cell])^2)
                derivative = @. profile * (velocity - vline[cell]) / sigma[cell]^2 *
                    C_LIGHT_KMS / rest_frequency_Hz
                Iprofile .+= profile
                Vprofile .+= 0.5splitting_coefficient .* Bline[cell] .* derivative
            end
            Ipeak_map[pixel] = maximum(Iprofile)
            pV_map[pixel] = maximum(abs, Vprofile) /
                max(Ipeak_map[pixel], eps(Float64))
        end
        NHI_map, pV_map, Ipeak_map
    end

    zeeman_NHI_map, zeeman_pV_map, zeeman_Ipeak_map_K = zeeman_fraction_maps(
        zeeman_nHI, zeeman_v_components[los_dim], T, zeeman_Blos_uG_cube,
        los_dim, sky_dims, dx_los_cm, Float64(zeeman_microturbulence_kms),
        zeeman_nu0_Hz, zeeman_z, Int(zeeman_channel_count))
    zeeman_pV_numerator = apply_observational_beam_2d(
        zeeman_pV_map .* zeeman_Ipeak_map_K, cube, sky_dims)
    zeeman_Ipeak_map_K = apply_observational_beam_2d(
        zeeman_Ipeak_map_K, cube, sky_dims)
    zeeman_NHI_map = apply_observational_beam_2d(zeeman_NHI_map, cube, sky_dims)
    zeeman_pV_map = zeeman_pV_numerator ./ max.(zeeman_Ipeak_map_K, eps(Float64))
    fig_zeeman_p_column = polarization_column_figure(
        zeeman_NHI_map, zeeman_pV_map,
        L"100p_V\;[\%]", MHD_COLORS[2])
    display_zeeman_p_column ? fig_zeeman_p_column : nothing
end

# ╔═╡ fd817f74-8bc0-4df7-85ad-45a95522f80a
begin
    zeeman_true_text = @sprintf("%.5g", zeeman_Btrue_uG)
    zeeman_fit_text = @sprintf("%.5g", zeeman_Bfit_uG)
    zeeman_split_text = @sprintf("%.5g", zeeman_z * zeeman_Btrue_uG)
    Markdown.parse("""
    ### Selected-sightline Zeeman result

    | Quantity | Recovered value |
    |:--|:--|
    | ``\\mathrm{H\\,I}``-weighted ``B_{\\mathrm{LOS}}`` | **$(zeeman_true_text)** ``\\mu\\mathrm{G}`` |
    | Stokes-``V`` derivative fit | **$(zeeman_fit_text)** ``\\mu\\mathrm{G}`` |
    | ``\\mathrm{H\\,I}``-weighted splitting | **$(zeeman_split_text)** ``\\mathrm{Hz}`` |
    """)
end

# ╔═╡ a0030003-6f8c-4d0c-9a10-000000000003
begin
    zeeman_structure_specs = [
        (data = zeeman_Bmap_uG, label = L"\langle B_{\mathrm{LOS}}\rangle_{\mathrm{HI}}\;[\mu\mathrm{G}]", color = MHD_COLORS[1], period = nothing),
        (data = zeeman_split_map_Hz, label = L"\Delta\nu_Z\;[\mathrm{Hz}]", color = MHD_COLORS[2], period = nothing),
        (data = zeeman_NHI_map, label = L"N_{\mathrm{HI}}\;[\mathrm{cm}^{-2}]", color = MHD_COLORS[3], period = nothing),
        (data = zeeman_Ipeak_map_K, label = L"\max_v I(v)\;[\mathrm{K}]", color = MHD_COLORS[4], period = nothing),
        (data = zeeman_pV_map, label = L"p_V", color = MHD_COLORS[5], period = nothing),
    ]
    fig_zeeman_structure = display_observational_structure_functions ?
        observational_structure_figure(zeeman_structure_specs, cube, sky_dims,
            observational_structure_order, observational_structure_samples;
            heading = "Zeeman observable structure functions") : Figure(size = (900, 120))
    display_observational_structure_functions ? fig_zeeman_structure : nothing
end

# ╔═╡ 62b61ef2-8e5d-4fe9-a435-e18fb5be9461
md"""
---

## 16. MOOSE Faraday post-processing

Synchrotron emission, Faraday rotation, instrumental filtering, and RM synthesis.

| MOOSE figure | Display |
|:--|:--:|
| Faraday and synchrotron maps | $(@bind display_moose PlutoUI.CheckBox(default = true)) |
| Faraday tomography | $(@bind display_moose_tomography PlutoUI.CheckBox(default = true)) |
| Polarization fraction versus $N_{\rm H}$ | $(@bind display_moose_p_column PlutoUI.CheckBox(default = true)) |

| MOOSE setting | Control |
|:--|:--|
| Electron prescription | $(@bind moose_electron_model PlutoUI.Select(["Two-phase ionization", "Constant ionization fraction"]; default = "Two-phase ionization")) |
| Constant $x_e$ | $(@bind moose_constant_xe PlutoUI.NumberField(default = 0.01)) |
| CNM $x_e$ | $(@bind moose_cnm_xe PlutoUI.NumberField(default = 1.0e-4)) |
| WNM $x_e$ | $(@bind moose_wnm_xe PlutoUI.NumberField(default = 1.0e-2)) |
| CNM/WNM transition temperature [$\mathrm{K}$] | $(@bind moose_transition_T PlutoUI.NumberField(default = 200.0)) |
| Cosmic-ray electron index $p$ | $(@bind moose_cr_index PlutoUI.NumberField(1.0:0.1:5.0; default = 3.0)) |
| Observing frequency [$\mathrm{MHz}$] | $(@bind moose_frequency_MHz PlutoUI.NumberField(50.0:1.0:2000.0; default = 150.0)) |
| Synchrotron normalization [$\mathrm{K}\,(\mu\mathrm{G})^{-(p+1)/2}\,\mathrm{pc}^{-1}$] | $(@bind moose_synchrotron_norm PlutoUI.NumberField(default = 1.0)) |
| Apply MOOSE interferometric filtering | $(@bind apply_moose_interferometer PlutoUI.CheckBox(default = false)) |
| Largest retained Fourier scale [pixels] | $(@bind moose_largest_scale_pix PlutoUI.NumberField(2.0:1.0:4096.0; default = 154.0)) |
| Smallest retained Fourier scale [pixels] | $(@bind moose_smallest_scale_pix PlutoUI.NumberField(2.0:0.5:256.0; default = 2.0)) |
| Add MOOSE Gaussian noise to $Q/U$ | $(@bind add_moose_noise PlutoUI.CheckBox(default = false)) |
| Polarized signal-to-noise ratio | $(@bind moose_instrument_snr PlutoUI.NumberField(0.1:0.1:1000.0; default = 10.0)) |
| Instrument random seed | $(@bind moose_instrument_seed PlutoUI.NumberField(0:1:100000; default = 42)) |
| Display shifted $uv$ transfer mask | $(@bind show_moose_uv_mask PlutoUI.CheckBox(default = false)) |
| Faraday-depth map | $(@bind show_moose_phi PlutoUI.CheckBox(default = true)) |
| Synchrotron brightness | $(@bind show_moose_I PlutoUI.CheckBox(default = true)) |
| Polarized brightness | $(@bind show_moose_P PlutoUI.CheckBox(default = true)) |
| Polarization fraction | $(@bind show_moose_fraction PlutoUI.CheckBox(default = true)) |
| RM-synthesis band start [$\mathrm{MHz}$] | $(@bind moose_band_start_MHz PlutoUI.NumberField(30.0:1.0:2000.0; default = 120.0)) |
| RM-synthesis band end [$\mathrm{MHz}$] | $(@bind moose_band_end_MHz PlutoUI.NumberField(30.0:1.0:2000.0; default = 167.0)) |
| Frequency-channel width [$\mathrm{MHz}$] | $(@bind moose_band_step_MHz PlutoUI.NumberField(0.05:0.05:20.0; default = 1.0)) |
| Minimum Faraday depth [$\mathrm{rad\,m^{-2}}$] | $(@bind moose_phi_min PlutoUI.NumberField(-500.0:0.25:0.0; default = -20.0)) |
| Maximum Faraday depth [$\mathrm{rad\,m^{-2}}$] | $(@bind moose_phi_max PlutoUI.NumberField(0.0:0.25:500.0; default = 20.0)) |
| Faraday-depth step [$\mathrm{rad\,m^{-2}}$] | $(@bind moose_dphi PlutoUI.NumberField(0.05:0.05:10.0; default = 0.25)) |
| First sky-axis pixel for $F(\phi)$ | $(@bind moose_sky_i PlutoUI.Slider(1:size(cube.rho, sky_dims[1]); default = cld(size(cube.rho, sky_dims[1]), 2), show_value = true)) |
| Second sky-axis pixel for $F(\phi)$ | $(@bind moose_sky_j PlutoUI.Slider(1:size(cube.rho, sky_dims[2]); default = cld(size(cube.rho, sky_dims[2]), 2), show_value = true)) |
| Peak Faraday-spectrum map (pmax) | $(@bind show_moose_pmax PlutoUI.CheckBox(default = true)) |
| Faraday-spectrum amplitude | $(@bind show_moose_F_abs PlutoUI.CheckBox(default = true)) |
| $\Re F(\phi)$ spectrum | $(@bind show_moose_F_real PlutoUI.CheckBox(default = true)) |
| $\Im F(\phi)$ spectrum | $(@bind show_moose_F_imag PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 0e8d9cab-aef2-42cd-959d-973764340f08
begin
    moose_ne = if moose_electron_model == "Constant ionization fraction"
        Float64(moose_constant_xe) .* number_density_cells
    else
        ionization_fraction = ifelse.(T .< Float64(moose_transition_T),
            Float64(moose_cnm_xe), Float64(moose_wnm_xe))
        ionization_fraction .* number_density_cells
    end

    moose_Blos_uG = GAUSS_TO_MICROGAUSS .* Bcomponents[los_dim]
    moose_Bsky1_uG = GAUSS_TO_MICROGAUSS .* Bcomponents[sky_dims[1]]
    moose_Bsky2_uG = GAUSS_TO_MICROGAUSS .* Bcomponents[sky_dims[2]]
    moose_Bperp_uG = sqrt.(moose_Bsky1_uG .^ 2 .+ moose_Bsky2_uG .^ 2)
    moose_phi_increment = 0.812 .* moose_ne .* moose_Blos_uG .* dx_los_pc
    moose_phi_to_cell = cumsum(moose_phi_increment; dims = los_dim) .- 0.5 .* moose_phi_increment
    moose_phi_map = finite_sum_dims(moose_phi_increment, los_dim)

    moose_p = Float64(moose_cr_index)
    moose_B_exponent = (moose_p + 1) / 2
    moose_temperature_spectral_index = -(moose_p + 3) / 2
    moose_frequency_scale = (Float64(moose_frequency_MHz) / 150.0)^moose_temperature_spectral_index
    moose_emissivity_Kpc = Float64(moose_synchrotron_norm) .* moose_frequency_scale .*
        moose_Bperp_uG .^ moose_B_exponent
    moose_intrinsic_angle = atan.(moose_Bsky2_uG, moose_Bsky1_uG) .+ pi / 2
    moose_lambda2_m2 = (299_792_458.0 / (Float64(moose_frequency_MHz) * 1.0e6))^2
    moose_polarization_phase = 2 .* (moose_intrinsic_angle .+ moose_phi_to_cell .* moose_lambda2_m2)
    moose_I_K = finite_sum_dims(moose_emissivity_Kpc, los_dim) .* dx_los_pc
    moose_Q_K = finite_sum_dims(moose_emissivity_Kpc .* cos.(moose_polarization_phase),
        los_dim) .* dx_los_pc
    moose_U_K = finite_sum_dims(moose_emissivity_Kpc .* sin.(moose_polarization_phase),
        los_dim) .* dx_los_pc
    moose_uv_transfer = moose_instrument_transfer(size(moose_I_K),
        moose_largest_scale_pix, moose_smallest_scale_pix)
    if apply_moose_interferometer
        moose_I_K = apply_moose_interferometer_2d(moose_I_K, moose_uv_transfer)
        moose_Q_K = apply_moose_interferometer_2d(moose_Q_K, moose_uv_transfer)
        moose_U_K = apply_moose_interferometer_2d(moose_U_K, moose_uv_transfer)
    end
    moose_I_K = apply_observational_beam_2d(moose_I_K, cube, sky_dims)
    moose_Q_K = apply_observational_beam_2d(moose_Q_K, cube, sky_dims)
    moose_U_K = apply_observational_beam_2d(moose_U_K, cube, sky_dims)
    add_moose_noise && add_moose_qu_noise!(moose_Q_K, moose_U_K,
        moose_instrument_snr, MersenneTwister(Int(moose_instrument_seed)))
    moose_P_K = sqrt.(moose_Q_K .^ 2 .+ moose_U_K .^ 2)
    moose_fraction = moose_P_K ./ max.(abs.(moose_I_K), eps(Float64))

    moose_specs = NamedTuple[]
    show_moose_phi && push!(moose_specs, (data = moose_phi_map,
        label = L"\phi\;[\mathrm{rad\,m}^{-2}]", colormap = :balance, diverging = true))
    show_moose_I && push!(moose_specs, (data = moose_I_K,
        label = L"T_{\mathrm{syn}}\;[\mathrm{K}]", colormap = :magma, diverging = false))
    show_moose_P && push!(moose_specs, (data = moose_P_K,
        label = L"P_\nu\;[\mathrm{K}]", colormap = :viridis, diverging = false))
    show_moose_fraction && push!(moose_specs, (data = moose_fraction,
        label = L"P_\nu/I_\nu", colormap = :plasma, diverging = false))
    show_moose_uv_mask && push!(moose_specs, (data = FFTW.fftshift(moose_uv_transfer),
        label = L"H(u,v)", colormap = :grays, diverging = false))

    if isempty(moose_specs)
        fig_moose = Figure(size = (900, 180))
        Label(fig_moose[1, 1], L"\mathrm{Select\ at\ least\ one\ MOOSE\ product.}", fontsize = 20)
    else
        moose_ncols = length(moose_specs) == 1 ? 1 : 2
        moose_nrows = cld(length(moose_specs), moose_ncols)
        fig_moose = Figure(size = (560moose_ncols, 420moose_nrows))
        for (index, spec) in enumerate(moose_specs)
            row, col = cld(index, moose_ncols), mod1(index, moose_ncols)
            panel = fig_moose[row, col] = GridLayout()
            ax = latex_axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            moose_colorrange = robust_colorrange(spec.data, color_percentile; diverging = spec.diverging)
            hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap, colorrange = moose_colorrange)
            latex_colorbar(panel[1, 2], hm; label = as_latex(spec.label), tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_moose ? fig_moose : nothing
end

# ╔═╡ 9a0f6ce3-9d80-46ca-ba3d-be5c31eaf032
begin
    fig_moose_p_column = polarization_column_figure(
        column_density, moose_fraction,
        L"100p_{\mathrm F}=100P_\nu/I_\nu\;[\%]", MHD_COLORS[5])
    display_moose_p_column ? fig_moose_p_column : nothing
end

# ╔═╡ c734b8e0-0bf7-42fc-bd26-0a451dd5f5f7
begin
    moose_band_lo = min(Float64(moose_band_start_MHz), Float64(moose_band_end_MHz))
    moose_band_hi = max(Float64(moose_band_start_MHz), Float64(moose_band_end_MHz))
    moose_band_step = max(Float64(moose_band_step_MHz), 0.05)
    moose_band_frequency_MHz = collect(moose_band_lo:moose_band_step:moose_band_hi)
    length(moose_band_frequency_MHz) >= 2 ||
        (moose_band_frequency_MHz = [moose_band_lo, moose_band_hi + moose_band_step])
    moose_band_frequency_Hz = moose_band_frequency_MHz .* 1.0e6
    moose_band_lambda2_m2 = (299_792_458.0 ./ moose_band_frequency_Hz) .^ 2
    moose_phi_lo = min(Float64(moose_phi_min), Float64(moose_phi_max))
    moose_phi_hi = max(Float64(moose_phi_min), Float64(moose_phi_max))
    moose_phi_step = max(Float64(moose_dphi), 0.05)
    moose_phi_axis = collect(moose_phi_lo:moose_phi_step:moose_phi_hi)
    length(moose_phi_axis) >= 2 || (moose_phi_axis = [moose_phi_lo, moose_phi_lo + moose_phi_step])

    moose_sky_shape = size(moose_phi_map)
    moose_nfrequency = length(moose_band_frequency_MHz)
    moose_Q_band_K = Array{Float64}(undef, moose_sky_shape..., moose_nfrequency)
    moose_U_band_K = similar(moose_Q_band_K)
    moose_emissivity_base_Kpc = Float64(moose_synchrotron_norm) .* moose_Bperp_uG .^ moose_B_exponent
    for channel in eachindex(moose_band_frequency_MHz)
        frequency_scale = (moose_band_frequency_MHz[channel] / 150.0)^moose_temperature_spectral_index
        phase = 2 .* (moose_intrinsic_angle .+
            moose_phi_to_cell .* moose_band_lambda2_m2[channel])
        moose_Q_band_K[:, :, channel] .= finite_sum_dims(
            moose_emissivity_base_Kpc .* frequency_scale .* cos.(phase), los_dim) .* dx_los_pc
        moose_U_band_K[:, :, channel] .= finite_sum_dims(
            moose_emissivity_base_Kpc .* frequency_scale .* sin.(phase), los_dim) .* dx_los_pc
    end
    if apply_moose_interferometer
        moose_Q_band_K = apply_moose_interferometer_cube(
            moose_Q_band_K, moose_uv_transfer)
        moose_U_band_K = apply_moose_interferometer_cube(
            moose_U_band_K, moose_uv_transfer)
    end
    moose_Q_band_K = apply_observational_beam_cube(
        moose_Q_band_K, cube, sky_dims)
    moose_U_band_K = apply_observational_beam_cube(
        moose_U_band_K, cube, sky_dims)
    add_moose_noise && add_moose_qu_noise!(moose_Q_band_K, moose_U_band_K,
        moose_instrument_snr, MersenneTwister(Int(moose_instrument_seed)))

    moose_lambda0_sq_m2 = mean(moose_band_lambda2_m2)
    moose_rm_phase_matrix = [cis(-2.0 * phi * (lambda2 - moose_lambda0_sq_m2))
        for lambda2 in moose_band_lambda2_m2, phi in moose_phi_axis]
    moose_P_band_matrix = reshape(complex.(moose_Q_band_K, moose_U_band_K), :, moose_nfrequency)
    moose_F_matrix = moose_P_band_matrix * moose_rm_phase_matrix / moose_nfrequency
    moose_F_complex = reshape(moose_F_matrix, moose_sky_shape..., length(moose_phi_axis))
    moose_F_abs = abs.(moose_F_complex)
    moose_pmax_K = finite_maximum_dims(moose_F_abs, 3)
end

# ╔═╡ e9c46999-b6cb-4bf4-93ff-e23d727698e1
begin
    moose_tomography_specs = Symbol[]
    show_moose_pmax && push!(moose_tomography_specs, :pmax)
    (show_moose_F_abs || show_moose_F_real || show_moose_F_imag) &&
        push!(moose_tomography_specs, :spectrum)
    if isempty(moose_tomography_specs)
        fig_moose_tomography = Figure(size = (900, 180))
        Label(fig_moose_tomography[1, 1],
            L"\mathrm{Select\ the\ }p_{\max}\mathrm{\ map\ or\ an\ }F(\phi)\mathrm{\ component.}", fontsize = 20)
    else
        fig_moose_tomography = Figure(size = (570length(moose_tomography_specs), 410))
        for (index, product) in enumerate(moose_tomography_specs)
            if product == :pmax
                panel = fig_moose_tomography[1, index] = GridLayout()
                ax = latex_axis(panel[1, 1], xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                    ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
                hm = heatmap!(ax, sky_coordinates[1], sky_coordinates[2], moose_pmax_K;
                    colormap = :viridis,
                    colorrange = robust_colorrange(moose_pmax_K, color_percentile))
                scatter!(ax, [sky_coordinates[1][Int(moose_sky_i)]],
                    [sky_coordinates[2][Int(moose_sky_j)]];
                    marker = :cross, markersize = 20, strokewidth = 3, color = :white)
                latex_colorbar(panel[1, 2], hm; label = L"p_{\max}=\max_\phi|F(\phi)|\;[\mathrm{K}]",
                    tickformat = latex_ticklabels)
                colsize!(panel, 2, 22)
            else
                ax = latex_axis(fig_moose_tomography[1, index],
                    xlabel = L"\phi\;[\mathrm{rad\,m}^{-2}]", ylabel = L"F(\phi)\;[\mathrm{K}]")
                spectrum = @view moose_F_complex[Int(moose_sky_i), Int(moose_sky_j), :]
                show_moose_F_abs && lines!(ax, moose_phi_axis, abs.(spectrum);
                    color = :black, linewidth = 2.8, label = "|F(φ)|")
                show_moose_F_real && lines!(ax, moose_phi_axis, real.(spectrum);
                    color = MHD_COLORS[1], linewidth = 2, label = "Re F(φ)")
                show_moose_F_imag && lines!(ax, moose_phi_axis, imag.(spectrum);
                    color = MHD_COLORS[2], linewidth = 2, label = "Im F(φ)")
                hlines!(ax, [0.0]; color = (:gray, 0.5), linestyle = :dot)
                axislegend(ax; position = :rt, framevisible = false)
            end
        end
    end
    display_moose_tomography ? fig_moose_tomography : nothing
end

# ╔═╡ a0040004-6f8c-4d0c-9a10-000000000004
begin
    moose_structure_specs = [
        (data = moose_phi_map, label = L"\phi\;[\mathrm{rad\,m}^{-2}]", color = MHD_COLORS[1], period = nothing),
        (data = moose_I_K, label = L"T_{\mathrm{syn}}\;[\mathrm{K}]", color = MHD_COLORS[2], period = nothing),
        (data = moose_Q_K, label = L"Q_\nu\;[\mathrm{K}]", color = MHD_COLORS[3], period = nothing),
        (data = moose_U_K, label = L"U_\nu\;[\mathrm{K}]", color = MHD_COLORS[4], period = nothing),
        (data = moose_P_K, label = L"P_\nu\;[\mathrm{K}]", color = MHD_COLORS[5], period = nothing),
        (data = moose_fraction, label = L"P_\nu/I_\nu", color = MHD_COLORS[6], period = nothing),
        (data = moose_pmax_K, label = L"p_{\max}\;[\mathrm{K}]", color = MHD_COLORS[1], period = nothing),
    ]
    fig_moose_structure = display_observational_structure_functions ?
        observational_structure_figure(moose_structure_specs, cube, sky_dims,
            observational_structure_order, observational_structure_samples;
            heading = "MOOSE observable structure functions") : Figure(size = (900, 120))
    display_observational_structure_functions ? fig_moose_structure : nothing
end

# ╔═╡ ab1a7df4-ae91-47db-cb4e-cf6d42fb0143
md"""
---

## 17. Polarization fraction versus intensity

These two-dimensional histograms show the sky-pixel distribution of polarization fraction versus total intensity. Empty bins are masked and populated bins are colored by their probability mass. The displayed fraction is $P/I$ for dust, dichroic starlight, and Faraday synchrotron emission. For Zeeman splitting, the analogous circular fraction is $p_V=\max_v|V(v)|/\max_v I(v)$.

| Two-dimensional histogram | Display |
|:--|:--:|
| Display selected histograms | $(@bind display_polarization_intensity_histograms PlutoUI.CheckBox(default = true)) |
| Dust $P_\nu/I_\nu$ | $(@bind show_dust_p_intensity_histogram PlutoUI.CheckBox(default = true)) |
| Starlight $P_\star/I_\star$ | $(@bind show_starlight_p_intensity_histogram PlutoUI.CheckBox(default = false)) |
| Faraday $P_\nu/I_\nu$ | $(@bind show_faraday_p_intensity_histogram PlutoUI.CheckBox(default = false)) |
| Zeeman $p_V$ | $(@bind show_zeeman_p_intensity_histogram PlutoUI.CheckBox(default = false)) |

| Histogram setting | Control |
|:--|:--|
| Number of bins per axis | $(@bind polarization_intensity_bins PlutoUI.Slider(20:5:100; default = 50, show_value = true)) |
| Logarithmic intensity axis | $(@bind log_polarization_intensity PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ bc2b8e05-bfa2-48ec-dc5f-d07e530c1254
begin
    function polarization_intensity_histogram(intensity, fraction, bin_count, use_log)
        I = vec(Float64.(intensity))
        p = 100 .* vec(Float64.(fraction))
        valid = isfinite.(I) .& isfinite.(p) .& (p .>= 0) .&
            (use_log ? (I .> 0) : trues(length(I)))
        I, p = I[valid], p[valid]
        isempty(I) && return (Float64[], Float64[], zeros(0, 0))
        I_low, I_high = length(I) >= 20 ? quantile(I, (0.001, 0.999)) : extrema(I)
        p_low, p_high = 0.0, length(p) >= 20 ? quantile(p, 0.999) : maximum(p)
        I_high > I_low || (I_high = I_low + max(abs(I_low), 1.0) * 1.0e-6)
        p_high > p_low || (p_high = p_low + max(abs(p_low), 1.0) * 1.0e-6)
        xedges = use_log ?
            10.0 .^ range(log10(max(I_low, floatmin(Float64))), log10(I_high);
                length = bin_count + 1) :
            collect(range(I_low, I_high; length = bin_count + 1))
        yedges = collect(range(p_low, p_high; length = bin_count + 1))
        I_clipped = clamp.(I, first(xedges), prevfloat(last(xedges)))
        p_clipped = clamp.(p, first(yedges), prevfloat(last(yedges)))
        histogram = fit(Histogram, (I_clipped, p_clipped), (xedges, yedges))
        total = sum(histogram.weights)
        log_probability = map(histogram.weights) do weight
            weight > 0 && total > 0 ? log10(weight / total) : NaN
        end
        xedges, yedges, log_probability
    end

    polarization_intensity_specs = NamedTuple[]
    show_dust_p_intensity_histogram && push!(polarization_intensity_specs, (
        intensity = dust_I, fraction = dust_fraction,
        xlabel = L"I_\nu^{\mathrm d}\;[\mathrm{MJy\,sr}^{-1}]",
        ylabel = L"100P_\nu^{\mathrm d}/I_\nu^{\mathrm d}\;[\%]"))
    show_starlight_p_intensity_histogram && push!(polarization_intensity_specs, (
        intensity = starlight_I_normalized, fraction = starlight_p_map,
        xlabel = L"I_\star/I_0",
        ylabel = L"100P_\star/I_\star\;[\%]"))
    show_faraday_p_intensity_histogram && push!(polarization_intensity_specs, (
        intensity = moose_I_K, fraction = moose_fraction,
        xlabel = L"I_\nu^{\mathrm F}\;[\mathrm{K}]",
        ylabel = L"100P_\nu^{\mathrm F}/I_\nu^{\mathrm F}\;[\%]"))
    show_zeeman_p_intensity_histogram && push!(polarization_intensity_specs, (
        intensity = zeeman_Ipeak_map_K, fraction = zeeman_pV_map,
        xlabel = L"\max_v I(v)\;[\mathrm{K}]",
        ylabel = L"100\max_v|V(v)|/\max_v I(v)\;[\%]"))

    if isempty(polarization_intensity_specs)
        fig_polarization_intensity = Figure(size = (900, 180))
        Label(fig_polarization_intensity[1, 1],
            L"\mathrm{Select\ at\ least\ one\ polarization\ histogram.}"; fontsize = 20)
    else
        panel_count = length(polarization_intensity_specs)
        panel_columns = panel_count == 1 ? 1 : 2
        panel_rows = cld(panel_count, panel_columns)
        fig_polarization_intensity = Figure(
            size = (570panel_columns, 430panel_rows))
        for (index, spec) in enumerate(polarization_intensity_specs)
            row, column = cld(index, panel_columns), mod1(index, panel_columns)
            panel = fig_polarization_intensity[row, column] = GridLayout()
            xedges, yedges, log_probability = polarization_intensity_histogram(
                spec.intensity, spec.fraction, Int(polarization_intensity_bins),
                log_polarization_intensity)
            if isempty(xedges)
                Label(panel[1, 1], L"\mathrm{No\ valid\ samples.}"; fontsize = 18)
                continue
            end
            axis = latex_axis(panel[1, 1], xlabel = spec.xlabel, ylabel = spec.ylabel,
                xscale = log_polarization_intensity ? log10 : identity,
                xticks = log_polarization_intensity ? DECADE_TICKS : Makie.automatic,
                xminorticks = log_polarization_intensity ? IntervalsBetween(9) : IntervalsBetween(5),
                xminorticksvisible = true)
            heat = heatmap!(axis, xedges, yedges, log_probability;
                colormap = :magma,
                colorrange = robust_colorrange(log_probability, color_percentile))
            latex_colorbar(panel[1, 2], heat;
                label = L"\log_{10}\mathcal{P}_{\mathrm{bin}}",
                tickformat = latex_ticklabels)
            colsize!(panel, 2, 22)
        end
    end
    display_polarization_intensity_histograms ? fig_polarization_intensity : nothing
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

# ╔═╡ 14e7606a-3a13-4c8e-b860-e40dc63a6fa2
md"""
---

## 19. Polarization fractions versus time

Each panel follows one synthetic-observation polarization fraction through the selected simulations. The reported value is either the area-weighted sky mean,

```math
\langle p\rangle_A=\frac{1}{N_{\rm pix}}\sum_j p_j,
```

or the intensity-weighted mean,

```math
\langle p\rangle_I=\frac{\sum_j P_j}{\sum_j I_j}=\frac{\sum_j p_jI_j}{\sum_j I_j}.
```

The thick horizontal segments above every panel identify two time intervals for the active run: the exponential dynamo interval used to fit $\Gamma_B$ in section 5 and the user-selected statistically steady interval. Thin vertical guides mark the end of growth and the beginning of the steady regime.

**Display polarization fractions versus time:** $(@bind display_polarization_time PlutoUI.CheckBox(default = false))

| Temporal polarization diagnostic | Display |
|:--|:--:|
| Thermal dust | $(@bind show_polarization_time_dust PlutoUI.CheckBox(default = true)) |
| Dichroic starlight | $(@bind show_polarization_time_starlight PlutoUI.CheckBox(default = true)) |
| Faraday synchrotron | $(@bind show_polarization_time_faraday PlutoUI.CheckBox(default = true)) |
| Zeeman circular polarization | $(@bind show_polarization_time_zeeman PlutoUI.CheckBox(default = true)) |

| Temporal setting | Control |
|:--|:--|
| Sky averaging | $(@bind polarization_time_weighting PlutoUI.Select(["Intensity-weighted mean" => "intensity", "Area-weighted mean" => "area"]; default = "intensity")) |
| Snapshot stride | $(@bind polarization_time_stride PlutoUI.Slider(1:maximum_snapshot_count; default = 1, show_value = true)) |
| Zeeman velocity channels | $(@bind polarization_time_zeeman_channels PlutoUI.Select([31, 51, 81, 101]; default = 51)) |
| Display growth and steady-state bars | $(@bind show_polarization_regime_bars PlutoUI.CheckBox(default = true)) |
| Steady-state start snapshot | $(@bind polarization_steady_start PlutoUI.Slider(1:length(run_files[selected_run]); default = min(last(growth_fit_window) + 1, length(run_files[selected_run])), show_value = true)) |
"""

# ╔═╡ b37d82bc-7819-47f9-b0ab-4de2f124f3cc
begin
    fig_polarization_time = let
    function polarization_sky_average(fraction, intensity, weighting)
        p = vec(Float64.(fraction))
        I = vec(Float64.(intensity))
        valid = isfinite.(p) .& isfinite.(I) .& (p .>= 0) .& (I .>= 0)
        any(valid) || return NaN
        if weighting == "intensity"
            normalization = sum(I[valid])
            normalization > 0 || return NaN
            return sum(p[valid] .* I[valid]) / normalization
        end
        mean(p[valid])
    end

    function temporal_dust_fraction(c, line_dim, plane_dims, weighting)
        components = (c.bx, c.by, c.bz)
        B2 = c.bx .^ 2 .+ c.by .^ 2 .+ c.bz .^ 2
        B1, B2sky = components[plane_dims[1]], components[plane_dims[2]]
        cos2gamma = clamp.((B1 .^ 2 .+ B2sky .^ 2) ./ max.(B2, eps(Float64)), 0.0, 1.0)
        angle = atan.(B2sky, B1) .+ pi / 2
        dx_cm = c.L[line_dim] / size(c.rho, line_dim) * PC_CM
        nH = c.rho ./ (max(Float64(dust_mu_H), eps(Float64)) * M_H_CGS)
        column_weight = dx_cm .* nH
        p0 = Float64(dust_p0)
        I = finite_sum_dims(column_weight .* (1 .- p0 .* (cos2gamma .- 2 / 3)),
            line_dim)
        Q = finite_sum_dims(column_weight .* p0 .* cos2gamma .* cos.(2 .* angle),
            line_dim)
        U = finite_sum_dims(column_weight .* p0 .* cos2gamma .* sin.(2 .* angle),
            line_dim)
        I = apply_observational_beam_2d(I, c, plane_dims)
        Q = apply_observational_beam_2d(Q, c, plane_dims)
        U = apply_observational_beam_2d(U, c, plane_dims)
        P = sqrt.(Q .^ 2 .+ U .^ 2)
        polarization_sky_average(P ./ max.(I, eps(Float64)), I, weighting)
    end

    function temporal_starlight_fraction(c, line_dim, plane_dims, weighting)
        nlos = size(c.rho, line_dim)
        dx_pc = c.L[line_dim] / nlos
        cell_count = clamp(ceil(Int, clamp(Float64(starlight_star_distance_pc),
            dx_pc, c.L[line_dim]) / dx_pc), 1, nlos)
        map_shape = size(selectdim(c.rho, line_dim, 1))
        I = fill(Float64(starlight_initial_I), map_shape)
        Q = fill(Float64(starlight_initial_Q), map_shape)
        U = fill(Float64(starlight_initial_U), map_shape)
        V = fill(Float64(starlight_initial_V), map_shape)
        components = (c.bx, c.by, c.bz)
        p0 = clamp(Float64(starlight_p0), 0.0, 1 - eps(Float64))
        nh_per_av = max(Float64(starlight_nh_per_av), eps(Float64))
        for cell_index in 1:cell_count
            density = selectdim(c.rho, line_dim, cell_index) ./
                (max(Float64(starlight_mu_H), eps(Float64)) * M_H_CGS)
            east = selectdim(components[plane_dims[1]], line_dim, cell_index)
            north = selectdim(components[plane_dims[2]], line_dim, cell_index)
            psi = atan.(east, north)
            tau = density .* (dx_pc * PC_CM) ./ nh_per_av
            I, Q, U, V = starlight_mueller_step(I, Q, U, V, tau, psi, p0)
        end
        I = apply_observational_beam_2d(I, c, plane_dims)
        Q = apply_observational_beam_2d(Q, c, plane_dims)
        U = apply_observational_beam_2d(U, c, plane_dims)
        fraction = clamp.(sqrt.(Q .^ 2 .+ U .^ 2) ./
            max.(abs.(I), eps(Float64)), 0.0, 1.0)
        polarization_sky_average(fraction, abs.(I), weighting)
    end

    function temporal_faraday_fraction(c, line_dim, plane_dims, weighting)
        components = (c.bx, c.by, c.bz)
        B1 = GAUSS_TO_MICROGAUSS .* components[plane_dims[1]]
        B2sky = GAUSS_TO_MICROGAUSS .* components[plane_dims[2]]
        Blos = GAUSS_TO_MICROGAUSS .* components[line_dim]
        Bperp = sqrt.(B1 .^ 2 .+ B2sky .^ 2)
        temperature = Float64(mean_molecular_weight) * M_H_CGS .* c.P ./
            (K_B_CGS .* c.rho)
        number_density_local = c.rho ./ (Float64(mean_molecular_weight) * M_H_CGS)
        electron_density = if moose_electron_model == "Constant ionization fraction"
            Float64(moose_constant_xe) .* number_density_local
        else
            ionization = ifelse.(temperature .< Float64(moose_transition_T),
                Float64(moose_cnm_xe), Float64(moose_wnm_xe))
            ionization .* number_density_local
        end
        dx_pc = c.L[line_dim] / size(c.rho, line_dim)
        phi_increment = 0.812 .* electron_density .* Blos .* dx_pc
        phi_to_cell = cumsum(phi_increment; dims = line_dim) .- 0.5 .* phi_increment
        cr_index = Float64(moose_cr_index)
        exponent = (cr_index + 1) / 2
        spectral_index = -(cr_index + 3) / 2
        frequency_scale = (Float64(moose_frequency_MHz) / 150.0)^spectral_index
        emissivity = Float64(moose_synchrotron_norm) .* frequency_scale .* Bperp .^ exponent
        intrinsic_angle = atan.(B2sky, B1) .+ pi / 2
        lambda2 = (299_792_458.0 / (Float64(moose_frequency_MHz) * 1.0e6))^2
        phase = 2 .* (intrinsic_angle .+ phi_to_cell .* lambda2)
        I = finite_sum_dims(emissivity, line_dim) .* dx_pc
        Q = finite_sum_dims(emissivity .* cos.(phase), line_dim) .* dx_pc
        U = finite_sum_dims(emissivity .* sin.(phase), line_dim) .* dx_pc
        transfer = moose_instrument_transfer(size(I), moose_largest_scale_pix,
            moose_smallest_scale_pix)
        if apply_moose_interferometer
            I = apply_moose_interferometer_2d(I, transfer)
            Q = apply_moose_interferometer_2d(Q, transfer)
            U = apply_moose_interferometer_2d(U, transfer)
        end
        I = apply_observational_beam_2d(I, c, plane_dims)
        Q = apply_observational_beam_2d(Q, c, plane_dims)
        U = apply_observational_beam_2d(U, c, plane_dims)
        if add_moose_noise
            temporal_seed = Int(moose_instrument_seed) +
                round(Int, 1000abs(Float64(c.t)))
            add_moose_qu_noise!(Q, U, moose_instrument_snr,
                MersenneTwister(temporal_seed))
        end
        P = sqrt.(Q .^ 2 .+ U .^ 2)
        polarization_sky_average(P ./ max.(abs.(I), eps(Float64)), abs.(I), weighting)
    end

    function temporal_zeeman_fraction(c, line_dim, plane_dims, weighting, channel_count)
        components = (c.vx, c.vy, c.vz)
        magnetic_components = (c.bx, c.by, c.bz)
        temperature = Float64(mean_molecular_weight) * M_H_CGS .* c.P ./
            (K_B_CGS .* c.rho)
        nHI = Float64(zeeman_neutral_fraction) .* c.rho ./ (1.4M_H_CGS)
        Blos = GAUSS_TO_MICROGAUSS .* magnetic_components[line_dim]
        dx_cm = c.L[line_dim] / size(c.rho, line_dim) * PC_CM
        _, fraction, Ipeak = zeeman_fraction_maps(nHI, components[line_dim],
            temperature, Blos, line_dim, plane_dims, dx_cm,
            Float64(zeeman_microturbulence_kms),
            Float64(zeeman_frequency_MHz) * 1.0e6,
            Float64(zeeman_coefficient_Hz_uG), channel_count)
        polarized_peak = apply_observational_beam_2d(fraction .* Ipeak, c, plane_dims)
        Ipeak = apply_observational_beam_2d(Ipeak, c, plane_dims)
        fraction = polarized_peak ./ max.(Ipeak, eps(Float64))
        polarization_sky_average(fraction, Ipeak, weighting)
    end

    polarization_time_panel_specs = NamedTuple[]
    show_polarization_time_dust && push!(polarization_time_panel_specs,
        (field = :dust, ylabel = L"100\langle p_{\mathrm d}\rangle\;[\%]"))
    show_polarization_time_starlight && push!(polarization_time_panel_specs,
        (field = :starlight, ylabel = L"100\langle p_{\star}\rangle\;[\%]"))
    show_polarization_time_faraday && push!(polarization_time_panel_specs,
        (field = :faraday, ylabel = L"100\langle p_{\mathrm F}\rangle\;[\%]"))
    show_polarization_time_zeeman && push!(polarization_time_panel_specs,
        (field = :zeeman, ylabel = L"100\langle p_V\rangle\;[\%]"))

    polarization_time_series = Dict{String, Any}()
    if display_polarization_time && !isempty(polarization_time_panel_specs)
        stride = max(Int(polarization_time_stride), 1)
        requested_fields = Set(getfield.(polarization_time_panel_specs, :field))
        for label in comparison_run_labels
            paths = run_files[label]
            indices = unique(vcat(collect(1:stride:length(paths)), length(paths)))
            values = NamedTuple[]
            for snapshot_index in indices
                local_cube = load_cube(paths[snapshot_index])
                push!(values, (
                    t = Float64(local_cube.t),
                    dust = :dust in requested_fields ?
                        temporal_dust_fraction(local_cube, los_dim, sky_dims,
                            polarization_time_weighting) : NaN,
                    starlight = :starlight in requested_fields ?
                        temporal_starlight_fraction(local_cube, los_dim, sky_dims,
                            polarization_time_weighting) : NaN,
                    faraday = :faraday in requested_fields ?
                        temporal_faraday_fraction(local_cube, los_dim, sky_dims,
                            polarization_time_weighting) : NaN,
                    zeeman = :zeeman in requested_fields ?
                        temporal_zeeman_fraction(local_cube, los_dim, sky_dims,
                            polarization_time_weighting,
                            Int(polarization_time_zeeman_channels)) : NaN,
                ))
            end
            polarization_time_series[label] = values
        end
    end

    if isempty(polarization_time_panel_specs)
        fig_polarization_time = Figure(size = (900, 180))
        Label(fig_polarization_time[1, 1],
            L"\mathrm{Select\ at\ least\ one\ temporal\ polarization\ diagnostic.}";
            fontsize = 20)
    elseif !display_polarization_time
        fig_polarization_time = Figure(size = (900, 180))
        Label(fig_polarization_time[1, 1],
            L"\mathrm{Temporal\ polarization\ diagnostics\ are\ disabled.}";
            fontsize = 20)
    else
        panel_columns = length(polarization_time_panel_specs) == 1 ? 1 : 2
        panel_rows = cld(length(polarization_time_panel_specs), panel_columns)
        fig_polarization_time = Figure(size = (560panel_columns, 400panel_rows + 95))
        growth_start_time = growth_has_interval ? growth_times[growth_first_index] : NaN
        growth_end_time = growth_has_interval ? growth_times[growth_last_index] : NaN
        steady_index = clamp(Int(polarization_steady_start), 1, length(growth_times))
        steady_start_time = growth_times[steady_index]
        steady_end_time = last(growth_times)

        for (panel_index, spec) in enumerate(polarization_time_panel_specs)
            row, column = cld(panel_index, panel_columns), mod1(panel_index, panel_columns)
            axis = latex_axis(fig_polarization_time[row, column],
                xlabel = L"t\;[\mathrm{Myr}]", ylabel = spec.ylabel)
            all_percentages = Float64[]
            for label in comparison_run_labels
                series = polarization_time_series[label]
                times = Float64.(getfield.(series, :t))
                percentages = 100 .* Float64.(getfield.(series, spec.field))
                valid = isfinite.(times) .& isfinite.(percentages)
                append!(all_percentages, percentages[valid])
                lines!(axis, times[valid], percentages[valid];
                    color = run_colors[label], linewidth = 2.4)
                scatter!(axis, times[valid], percentages[valid];
                    color = run_colors[label], markersize = 7)
            end

            if show_polarization_regime_bars && !isempty(all_percentages)
                data_min, data_max = extrema(all_percentages)
                span = max(data_max - data_min, max(abs(data_max), 1.0) * 0.08)
                growth_bar_y = data_max + 0.12span
                steady_bar_y = data_max + 0.28span
                if growth_has_interval && isfinite(growth_start_time) && isfinite(growth_end_time)
                    lines!(axis, [growth_start_time, growth_end_time],
                        [growth_bar_y, growth_bar_y]; color = MHD_COLORS[2],
                        linewidth = 7, linecap = :round)
                    vlines!(axis, [growth_end_time]; color = (MHD_COLORS[2], 0.55),
                        linewidth = 1.4, linestyle = :dash)
                end
                if steady_end_time >= steady_start_time
                    lines!(axis, [steady_start_time, steady_end_time],
                        [steady_bar_y, steady_bar_y]; color = MHD_COLORS[3],
                        linewidth = 7, linecap = :round)
                    vlines!(axis, [steady_start_time]; color = (MHD_COLORS[3], 0.55),
                        linewidth = 1.4, linestyle = :dash)
                end
                ylims!(axis, data_min - 0.08span, data_max + 0.42span)
            end
        end

        polarization_legend_layout =
            GridLayout(fig_polarization_time[panel_rows + 1, 1:panel_columns])
        Legend(polarization_legend_layout[1, 1],
            [LineElement(color = run_colors[label], linewidth = 2.5)
                for label in comparison_run_labels],
            legend_run_label.(comparison_run_labels), "Simulation";
            orientation = :horizontal, tellheight = true, framevisible = false,
            labelsize = 11)
        if show_polarization_regime_bars
            Legend(polarization_legend_layout[1, 2], [
                LineElement(color = MHD_COLORS[2], linewidth = 7),
                LineElement(color = MHD_COLORS[3], linewidth = 7),
            ], ["ΓB growth", "Steady state"], "Time regime";
                orientation = :horizontal, tellheight = true, framevisible = false,
                labelsize = 11)
        end
    end
    fig_polarization_time
    end
    stable_pluto_figure(display_polarization_time, fig_polarization_time)
end

# ╔═╡ b47871f8-13c9-41b8-9f42-f74e90bac653
begin
    export_figure_options = [
        "heatmaps" => "Projected heatmaps",
        "pdfs" => "Probability density functions",
        "phase_diagram" => "Pressure-density phase diagram",
        "time_evolution" => "Global time evolution",
        "phase_magnetic_time" => "Magnetic field by thermal phase",
        "magnetic_fit" => "Magnetic exponential fit",
        "growth_rate_relations" => "Gamma versus time and Mach numbers",
        "normalized_magnetic_relations" => "Multi-snapshot ln(B/B0) relations",
        "normalized_magnetic_field" => "Normalized magnetic field",
        "magnetic_density" => "Magnetic field versus density",
        "hro" => "Histogram of relative orientations",
        "hog" => "Histogram of oriented gradients",
        "energy_ratios" => "Energy ratios by density",
        "energy_time" => "Energy ratios versus time",
        "vorticity" => "Vorticity map",
        "enstrophy_density" => "Enstrophy by density and family parameter",
        "power_spectra" => "Power spectra",
        "structure_functions" => "Structure functions",
        "dust_polarization" => "Dust-polarization maps",
        "dust_structure" => "Dust observable structure functions",
        "dust_pixel_spectrum" => "Dust I, Q, U spectrum at selected pixel",
        "dust_statistics" => "Dust-polarization statistics",
        "dust_p_column" => "Dust polarization fraction versus NH",
        "starlight_maps" => "Dichroic starlight-polarization maps",
        "starlight_structure" => "Starlight observable structure functions",
        "starlight_profiles" => "Dichroic starlight sight-line profiles",
        "starlight_p_column" => "Starlight polarization fraction versus NH",
        "zeeman_maps" => "Zeeman-splitting maps",
        "zeeman_structure" => "Zeeman observable structure functions",
        "zeeman_spectra" => "Zeeman Stokes spectra",
        "zeeman_p_column" => "Zeeman circular-polarization fraction versus NHI",
        "moose" => "MOOSE post-processing",
        "moose_structure" => "MOOSE observable structure functions",
        "moose_tomography" => "MOOSE F(phi) and pmax",
        "moose_p_column" => "Faraday polarization fraction versus NH",
        "polarization_intensity" => "Polarization fraction versus intensity 2D histograms",
        "shine" => "SHINE H I maps",
        "shine_structure" => "SHINE observable structure functions",
        "shine_rgb" => "SHINE H I velocity RGB composite",
        "shine_spectrum" => "SHINE H I spectrum",
        "polarization_time" => "Polarization fractions versus time",
    ]
    export_figure_registry = Dict(
        "heatmaps" => fig_maps,
        "pdfs" => fig_pdf,
        "phase_diagram" => fig_phase,
        "time_evolution" => fig_time,
        "phase_magnetic_time" => fig_phase_B_time,
        "magnetic_fit" => fig_growth,
        "growth_rate_relations" => fig_gamma_relations,
        "normalized_magnetic_relations" => fig_normalized_B_relations,
        "normalized_magnetic_field" => fig_logB,
        "magnetic_density" => fig_bn,
        "hro" => fig_hro,
        "hog" => fig_hog,
        "energy_ratios" => fig_energy,
        "energy_time" => fig_energy_time,
        "vorticity" => fig_vorticity,
        "enstrophy_density" => fig_enstrophy_density,
        "power_spectra" => fig_spectra,
        "structure_functions" => fig_structure,
        "dust_polarization" => fig_dust,
        "dust_structure" => fig_dust_structure,
        "dust_pixel_spectrum" => fig_dust_pixel_spectrum,
        "dust_statistics" => fig_dust_statistics,
        "dust_p_column" => fig_dust_p_column,
        "starlight_maps" => fig_starlight_maps,
        "starlight_structure" => fig_starlight_structure,
        "starlight_profiles" => fig_starlight_profiles,
        "starlight_p_column" => fig_starlight_p_column,
        "zeeman_maps" => fig_zeeman_maps,
        "zeeman_structure" => fig_zeeman_structure,
        "zeeman_spectra" => fig_zeeman_spectra,
        "zeeman_p_column" => fig_zeeman_p_column,
        "moose" => fig_moose,
        "moose_structure" => fig_moose_structure,
        "moose_tomography" => fig_moose_tomography,
        "moose_p_column" => fig_moose_p_column,
        "polarization_intensity" => fig_polarization_intensity,
        "shine" => fig_shine,
        "shine_structure" => fig_shine_structure,
        "shine_rgb" => fig_shine_rgb,
        "shine_spectrum" => fig_shine_spectrum,
        "polarization_time" => fig_polarization_time,
    )
    nothing
end

# ╔═╡ 3ca92f4b-1f96-4f3e-9536-ce17ae786cc6
md"""
---

## 20. Figure export

1. Select the figure currently shown in the notebook.
2. Choose **PNG** for a raster image or **PDF** for vector output.
3. Use the download button below the table.

The export reproduces the active run, snapshot, line of sight, panel checkboxes, and slider values. PDF is recommended for papers because text and line art remain vectorial; PNG is convenient for slides and quick previews.

| Export setting | Control |
|:--|:--|
| Figure | $(@bind export_figure_key PlutoUI.Select(export_figure_options; default = "energy_ratios")) |
| Format | $(@bind export_figure_format PlutoUI.Select(["PNG", "PDF"]; default = "PDF")) |
"""

# ╔═╡ 1c93c944-f221-4055-843d-7ecefc2be6ed
begin
    export_extension = lowercase(export_figure_format)
    export_mime = export_figure_format == "PNG" ? MIME"image/png"() : MIME"application/pdf"()
    export_buffer = IOBuffer()
    show(export_buffer, export_mime, export_figure_registry[export_figure_key])
    export_bytes = take!(export_buffer)
    export_run_slug = replace(lowercase(selected_run), r"[^a-z0-9]+" => "_")
    export_filename = "mhd_$(export_figure_key)_$(export_run_slug)_snapshot_$(lpad(selected_snapshot, 3, '0')).$(export_extension)"
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
# ╟─3ef88702-bef5-4eca-a151-df97aa7ec2c4
# ╟─bdc11245-d76c-45fe-b79d-7b64861f5f53
# ╟─7bd6f2c9-ae49-4636-a251-f526ab347125
# ╟─55cbdbf4-e0f2-431e-a736-09f41ab7ee75
# ╟─98360288-85ca-4551-bdde-c12c7a329302
# ╟─6a9cfec2-2c80-4d72-94c0-cdb47aa4f046
# ╟─18e0d2d1-56ca-46cc-a13e-77f142612b5a
# ╟─3df9ad9a-b865-41f5-8d7a-34a52d0292dd
# ╟─290ecb0a-880a-44bd-9ddc-e6ee12b41a06
# ╟─e734297f-506e-45e1-8cb7-b2ae671893eb
# ╠═7174c31f-f186-48f9-b66e-29cf9a1c1fe3
# ╟─353cc6fb-c801-448a-a2c5-23dfd1541704
# ╟─a8ef96ab-0ddd-4eb2-a216-b7d96c2a9a08
# ╟─32110739-b60e-4592-856a-dd74f7a37401
# ╟─5fbd0d53-09e8-4b39-bb1d-29c8cc45c6ee
# ╟─94a0a0dc-baf6-4e62-a51e-dc6124d98fd4
# ╟─c12d1f54-40b8-4865-9562-8dcb519f924a
# ╠═76d249f9-d6fd-4513-aa76-7fb386058c37
# ╟─72626861-42ce-4ac0-b980-78f498f8a629
# ╟─496cbf2d-77a1-4a1a-b760-d4f8ea2ea9de
# ╟─1e0e6c0e-1ae0-40a1-a6f8-9c18fae91961
# ╟─fa155a62-da75-4530-b5e9-215fd4f66412
# ╟─9da572fe-21f2-43df-9320-b8742fd64773
# ╠═41b4eb12-889d-43b3-87c2-fc7cccf8679f
# ╟─904ba663-d536-4b27-a379-4af927b0affb
# ╠═36aef377-3de7-435a-af83-3a90421e3159
# ╟─89c33295-34e7-49ec-8d04-52b2aac29cff
# ╟─8120f7e9-74ed-4a48-b2f4-dbc75ebf0132
# ╟─71c8ea26-d2ad-4430-9265-0b28d45bba1c
# ╠═f29d5dff-8810-4935-b880-9951bc1239fd
# ╟─cd7e037c-62f5-4682-92c2-92af7169d692
# ╠═7a4472c3-48c5-44c7-a487-d1021b84ee3a
# ╟─a7def784-6c52-4f98-9d76-ece082b45229
# ╟─dcc8f8f3-daaf-4d4b-92a3-4919ed5e36de
# ╠═3bed2efb-f849-4ddc-9f49-ed5e3d683370
# ╟─62bb58f1-37c3-4adb-a7fc-939f71a56635
# ╠═d56a8ba3-aa42-4351-a065-87275032a342
# ╟─1be45ca2-bd2f-4b77-9188-5b338d41483b
# ╟─89420b8d-c72e-4a04-91dd-043bc9ecef2e
# ╟─b1000001-6f8c-4d0c-9a10-000000000001
# ╠═b1000002-6f8c-4d0c-9a10-000000000002
# ╠═b1000003-6f8c-4d0c-9a10-000000000003
# ╠═b1000004-6f8c-4d0c-9a10-000000000004
# ╠═b1000005-6f8c-4d0c-9a10-000000000005
# ╟─298bd579-bb28-48b7-8c55-ea74804b9837
# ╟─e4a64ef5-9e6e-4ae2-8daf-42e4fc1724ba
# ╟─a1bf5e62-ecb8-47d7-8692-c33929c831ac
# ╟─5983ddd0-3ea9-4304-a1e2-390065309aef
# ╟─b7c04a6f-b24d-48ca-97d2-5a77d297aa8d
# ╟─2c4762f5-4780-4bfa-a320-e4102d9f5d6c
# ╟─72b2e359-c054-4efc-a834-a2f5fa99e3fc
# ╟─4e8ac6d4-4ee4-4e41-a63a-1c93213da33f
# ╟─d87379b7-3527-45a8-bc60-bec191c499af
# ╟─3190e127-1d53-49f1-bfab-b9645910c2c6
# ╟─df880704-b12e-49eb-84f2-67b6f9583a8a
# ╟─b02ceda6-0a0d-4543-91f8-a6669231ec72
# ╟─873f7ef2-719b-4ae6-b015-1a23c6c27836
# ╟─a8558c31-7dcf-433e-9950-a59e9acf158b
# ╟─3a731972-3404-478c-a572-00a05ab652b1
# ╠═24e60849-1c70-4df3-bd17-57d29949b7a6
# ╠═4b16d83f-4a1d-49e7-9270-8f573fd46835
# ╠═d69dd1ce-312a-48e0-a478-b470b299ed1b
# ╟─62440e86-b560-44ad-bb0a-43ae62e73fc3
# ╠═47b786d6-c7b5-44f4-946a-b8c485ad6380
# ╟─d6a2f4b1-59ac-4e77-a10a-4b74c0d89231
# ╟─130ccf03-d7cc-4b71-9210-dbc0e43dfa82
# ╟─5b3f6a91-246e-4cc3-8f68-164f7ff2f07c
# ╟─61c3b28c-9d62-4689-a85d-bc827e89641d
# ╠═a0010001-6f8c-4d0c-9a10-000000000001
# ╟─6f4e2d11-2a88-41f4-93dc-01b51d86fb4f
# ╠═a2d5319b-06f4-4efc-b3d7-3a9719292305
# ╟─d1e7626c-81d7-48cd-90be-695a13aa9997
# ╟─bcc9a131-9df2-455a-bdac-586323fd58d0
# ╟─0c8f7453-06e8-47ea-ae0f-4f09a89b16ae
# ╠═a0020002-6f8c-4d0c-9a10-000000000002
# ╟─7e8d4ac1-7b6e-44a8-981b-9c3a19c8de10
# ╟─67f95c39-1888-4d23-a2c2-2ee3a6cd7f0f
# ╟─76b9cf43-a13e-44c8-97d6-4760cc3aa486
# ╟─82e22c29-1cf5-47dc-a54c-68577f8069bc
# ╟─8f9e5bd2-8c7f-45b9-a92c-ad4b20d9ef21
# ╟─fd817f74-8bc0-4df7-85ad-45a95522f80a
# ╠═a0030003-6f8c-4d0c-9a10-000000000003
# ╟─62b61ef2-8e5d-4fe9-a435-e18fb5be9461
# ╟─0e8d9cab-aef2-42cd-959d-973764340f08
# ╟─9a0f6ce3-9d80-46ca-ba3d-be5c31eaf032
# ╟─c734b8e0-0bf7-42fc-bd26-0a451dd5f5f7
# ╟─e9c46999-b6cb-4bf4-93ff-e23d727698e1
# ╠═a0040004-6f8c-4d0c-9a10-000000000004
# ╟─ab1a7df4-ae91-47db-cb4e-cf6d42fb0143
# ╟─bc2b8e05-bfa2-48ec-dc5f-d07e530c1254
# ╟─478ec2f3-e057-4720-809c-17ca0a3dac21
# ╟─bcc05889-02bb-47cf-b672-139e8efe4137
# ╠═5ad762e1-105f-4cf7-9cf1-e0bb8c6f1bf5
# ╟─27e51ba5-4592-4766-9dde-0de383a889a0
# ╠═a0050005-6f8c-4d0c-9a10-000000000005
# ╟─14e7606a-3a13-4c8e-b860-e40dc63a6fa2
# ╟─b37d82bc-7819-47f9-b0ab-4de2f124f3cc
# ╠═b47871f8-13c9-41b8-9f42-f74e90bac653
# ╟─3ca92f4b-1f96-4f3e-9536-ce17ae786cc6
# ╠═1c93c944-f221-4055-843d-7ecefc2be6ed
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
