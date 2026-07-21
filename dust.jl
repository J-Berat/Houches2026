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

    const TICK_SUPERSCRIPTS = Dict(
        '-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³',
        '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
    )
    superscript_integer(value) =
        join(get(TICK_SUPERSCRIPTS, character, character) for character in string(value))

    """
    Render numeric ticks as one atomic TeX-styled text block.

    Genuine LaTeXStrings remain in axis and colorbar labels. Numeric tick arrays
    use Unicode mathematical glyphs because CairoMakie can otherwise update their
    LaTeX glyph blocks before the matching tick positions while Pluto is reactive.
    """
    function latex_number(x)
        if isnan(x)
            return "NaN"
        elseif isinf(x)
            return x > 0 ? "+∞" : "−∞"
        end
        x == 0 && return "0"
        exponent = floor(Int, log10(abs(x)))
        if exponent <= -3 || exponent >= 4
            mantissa = x / 10.0^exponent
            string(@sprintf("%.3g", mantissa), "×10", superscript_integer(exponent))
        else
            @sprintf("%.4g", x)
        end
    end

    latex_ticklabels(values) = latex_number.(values)
    as_latex(label::LaTeXString) = label
    as_latex(label::AbstractString) = latexstring(label)

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
# Dust — thermal polarization

Synthetic thermal-dust Stokes emission, polarization maps, spectra, statistics, and polarization fraction versus column density.

> **Reactive mode.** Open `dust.jl` with Pluto. Selecting a repository, run, snapshot, or line of sight updates all dependent products automatically.

> **Lazy startup.** Run `run_pluto.jl` with `DYNAMO_NOTEBOOK=dust.jl`. Pluto starts without evaluating the expensive cells; run the result cells you need and Pluto will resolve their upstream dependencies.

All dimensional quantities are converted to the physical units shown on their axes or colorbars. Projected means are density weighted unless stated otherwise, and periodic boundaries are used for spatial operations.
"""

# ╔═╡ bdc11245-d76c-45fe-b79d-7b64861f5f53
PlutoUI.TableOfContents(title = "Notebook sections", indent = true, depth = 3, aside = true)

# ╔═╡ 7bd6f2c9-ae49-4636-a251-f526ab347125
begin
    preferred_data_repository = get(
        ENV,
        "DYNAMO_DATA_REPOSITORY",
        "/Xnfs/Houches2026/DynSim/cooling_freq_output",
    )
    bundled_data_repository = joinpath(@__DIR__, "cooling", "VaryingMach")
    DEFAULT_DATA_REPOSITORY = isdir(preferred_data_repository) ?
        preferred_data_repository : bundled_data_repository
    SNAPSHOT_EXTENSIONS = (".h5", ".hdf5", ".fits", ".fit", ".fts")
    nothing
end

# ╔═╡ 55cbdbf4-e0f2-431e-a736-09f41ab7ee75
md"""
---

## 1. Data selection and physical units

### Data repository

1. Enter any repository, family, simulation, `DataCubes`, or direct snapshot-directory path.
2. Select **Load repository**. HDF5 and FITS runs and snapshots are discovered recursively.

| Data source | Control |
|:--|:--|
| Data path | $(@bind data_repository PlutoUI.confirm(PlutoUI.TextField(90; default = DEFAULT_DATA_REPOSITORY, placeholder = "/path/to/data"); label = "Load path")) |

Any folder name and nesting structure is accepted. FITS snapshots may be multi-extension files or directories containing one named FITS image per physical field.
"""

# ╔═╡ 98360288-85ca-4551-bdde-c12c7a329302
begin
    snapshot_extension(path) = lowercase(splitext(path)[2])
    is_snapshot_file(path) = isfile(path) && snapshot_extension(path) in SNAPSHOT_EXTENSIONS
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
            h5open(path, "r") do h
                any(dataset_path -> ndims(h[dataset_path]) == 3, hdf5_dataset_paths(h))
            end
        catch
            false
        end
    end

    function hdf5_field_path(h, field; required = true, source = "HDF5 file",
            overrides = Dict{Symbol,String}())
        aliases = Set(normalize_field_name.(FIELD_ALIASES[field]))
        paths = hdf5_dataset_paths(h)
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
            overrides = Dict{Symbol,String}())
        path = hdf5_field_path(h, field; required, source, overrides)
        isnothing(path) ? nothing : read(h[path])
    end

    hdf5_scalar_value(value) = value isa Number ? Float64(value) : Float64(only(value))

    function centered_hdf5_magnetic_component(h, centered, left, right, source;
            overrides = Dict{Symbol,String}())
        direct = read_hdf5_field(h, centered; required = false, source, overrides)
        isnothing(direct) || return direct
        lower = read_hdf5_field(h, left; required = false, source, overrides)
        upper = read_hdf5_field(h, right; required = false, source, overrides)
        (isnothing(lower) || isnothing(upper)) && error(
            "HDF5 magnetic component $(centered) in $(source) requires either " *
            "$(join(FIELD_ALIASES[centered], ", ")) or both face fields " *
            "$(join(FIELD_ALIASES[left], ", ")) and " *
            "$(join(FIELD_ALIASES[right], ", ")).",
        )
        0.5 .* (lower .+ upper)
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
            if hdf5_file_is_snapshot(path) || is_fits_file(path) ||
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
        !isempty(snapshot_sources(path)) && return [abspath(path)]
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
            "No HDF5 or FITS snapshots were found recursively under: $requested. " *
            "The folder may be a repository, family, run, DataCubes folder, or direct snapshot folder.",
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
                isnothing(hdu) ? nothing : Float64.(read(hdu))
            end
        else
            FITS(source) do fits
                hdu = fits_field_hdu(fits, field; primary_fallback)
                isnothing(hdu) ? nothing : Float64.(read(hdu))
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
        h5open(path, "r") do h
            candidates = vcat(names,
                ["metadata/$name" for name in names],
                ["parameters/$name" for name in names])
            for candidate in candidates
                haskey(h, candidate) || continue
                value = read(h[candidate])
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
        h5open(path, "r") do h
            density_path = hdf5_field_path(h, :rho; source = path,
                overrides)
            dimensions = size(h[density_path])
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
        h5open(hdf5_dataset_paths, HDF5_REFERENCE_FILE, "r")

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

        Reference file: `$(HDF5_REFERENCE_FILE)`

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

    function snapshot_files(label)
        sources = snapshot_sources(run_cube_directory(ROOT, RUN_DIRS[label]))
        isempty(sources) && error("No HDF5 or FITS snapshots found for $label in $ROOT")
        sources
    end

    function snapshot_time(path)
        if is_fits_file(path) || isdir(path)
            stored = read_fits_field(path, :t; required = false)
            !isnothing(stored) && length(stored) == 1 && return Float64(first(stored))
            header_time = fits_header_scalar(path, ["TIME", "T", "SIMTIME", "MYRTIME"])
            if isnothing(header_time)
                matched = match(r"([0-9]+(?:\.[0-9]+)?)\D*$", basename(path))
                return isnothing(matched) ? NaN : parse(Float64, matched.captures[1])
            end
            return header_time
        end
        h5open(path, "r") do h
            stored = read_hdf5_field(h, :t; required = false, source = path,
                overrides = HDF5_FIELD_OVERRIDES)
            if isnothing(stored)
                matched = match(r"([0-9]+(?:\.[0-9]+)?)\D*$", basename(path))
                return isnothing(matched) ? NaN : parse(Float64, matched.captures[1])
            end
            hdf5_scalar_value(stored)
        end
    end

    run_files = Dict(label => snapshot_files(label) for label in run_labels)
    run_times = Dict(label => snapshot_time.(run_files[label]) for label in run_labels)
    maximum_snapshot_count = maximum(length.(values(run_files)))
    run_summary = join([string(label, " = ", length(run_files[label]), " snapshots") for label in run_labels], ", ")
    nothing
end

# ╔═╡ e734297f-506e-45e1-8cb7-b2ae671893eb
md"""
| Navigation | Control |
|:--|:--|
| Run | $(@bind selected_run PlutoUI.Select(run_labels; default = run_labels[cld(length(run_labels), 2)])) |
| Maximum simulations opened | $(@bind maximum_open_cubes PlutoUI.NumberField(1:min(3, length(run_labels)); default = min(3, length(run_labels)))) |
| Simulations in comparative plots | $(@bind comparison_run_selection PlutoUI.MultiSelect(run_labels; default = run_labels[1:min(3, length(run_labels))])) |
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
        0.5 .* (lower .+ upper)
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

    "Load one HDF5 or FITS cube and center face-stored magnetic-field components."
    function load_cube(path)
        if is_fits_file(path) || isdir(path)
            raw_rho = read_fits_field(path, :rho; primary_fallback = true)
            raw_fields = (
                rho = raw_rho,
                P = read_fits_field(path, :P),
                vx = read_fits_field(path, :vx),
                vy = read_fits_field(path, :vy),
                vz = read_fits_field(path, :vz),
                bx = centered_fits_magnetic_component(path, :bx, :bx_l, :bx_r),
                by = centered_fits_magnetic_component(path, :by, :by_l, :by_r),
                bz = centered_fits_magnetic_component(path, :bz, :bz_l, :bz_r),
            )
            validate_cube_shapes(path, raw_fields)
            rho = Float64(density_unit_gcm3) .* raw_fields.rho
            P = Float64(pressure_unit_ergcm3) .* raw_fields.P
            velocity_scale = Float64(velocity_unit_kms)
            vx, vy, vz = velocity_scale .* raw_fields.vx,
                velocity_scale .* raw_fields.vy, velocity_scale .* raw_fields.vz
            magnetic_scale = Float64(magnetic_unit_G)
            bx, by, bz = magnetic_scale .* raw_fields.bx,
                magnetic_scale .* raw_fields.by, magnetic_scale .* raw_fields.bz
            L = Float64(length_unit_pc) .* fits_box_length(path, size(raw_rho))
            raw_time = snapshot_time(path)
            t = Float64(time_unit_Myr) * raw_time
            return (; rho, P, vx, vy, vz, bx, by, bz, L, t)
        end
        h5open(path, "r") do h
            raw_fields = (
                rho = read_hdf5_field(h, :rho; source = path,
                    overrides = HDF5_FIELD_OVERRIDES),
                P = read_hdf5_field(h, :P; source = path,
                    overrides = HDF5_FIELD_OVERRIDES),
                vx = read_hdf5_field(h, :vx; source = path,
                    overrides = HDF5_FIELD_OVERRIDES),
                vy = read_hdf5_field(h, :vy; source = path,
                    overrides = HDF5_FIELD_OVERRIDES),
                vz = read_hdf5_field(h, :vz; source = path,
                    overrides = HDF5_FIELD_OVERRIDES),
                bx = centered_hdf5_magnetic_component(h, :bx, :bx_l, :bx_r, path;
                    overrides = HDF5_FIELD_OVERRIDES),
                by = centered_hdf5_magnetic_component(h, :by, :by_l, :by_r, path;
                    overrides = HDF5_FIELD_OVERRIDES),
                bz = centered_hdf5_magnetic_component(h, :bz, :bz_l, :bz_r, path;
                    overrides = HDF5_FIELD_OVERRIDES),
            )
            validate_cube_shapes(path, raw_fields)
            rho = Float64(density_unit_gcm3) .* raw_fields.rho
            P = Float64(pressure_unit_ergcm3) .* raw_fields.P
            velocity_scale = Float64(velocity_unit_kms)
            vx, vy, vz = velocity_scale .* raw_fields.vx,
                velocity_scale .* raw_fields.vy, velocity_scale .* raw_fields.vz
            magnetic_scale = Float64(magnetic_unit_G)
            bx, by, bz = magnetic_scale .* raw_fields.bx,
                magnetic_scale .* raw_fields.by, magnetic_scale .* raw_fields.bz
            raw_length = read_hdf5_field(h, :L; source = path,
                overrides = HDF5_FIELD_OVERRIDES)
            length_values = raw_length isa Number ? fill(Float64(raw_length), 3) : Float64.(raw_length)
            length(length_values) == 3 || error(
                "HDF5 box length in $(path) must be a scalar or contain three values; " *
                "found size $(size(raw_length)).")
            L = Float64(length_unit_pc) .* length_values
            raw_time = read_hdf5_field(h, :t; required = false, source = path,
                overrides = HDF5_FIELD_OVERRIDES)
            stored_time = isnothing(raw_time) ? snapshot_time(path) : hdf5_scalar_value(raw_time)
            t = Float64(time_unit_Myr) .* stored_time
            (; rho, P, vx, vy, vz, bx, by, bz, L, t)
        end
    end

    number_density(rho) = rho ./ (Float64(mean_molecular_weight) * M_H_CGS)

    finite_values(A) = filter(isfinite, vec(Float64.(A)))
    finite_positive_values(A) = filter(x -> isfinite(x) && x > 0, vec(Float64.(A)))
    finite_mean(A; default = NaN) = begin
        values = finite_values(A)
        isempty(values) ? default : mean(values)
    end
    finite_quantile(A, q; default = NaN) = begin
        values = finite_values(A)
        isempty(values) ? default : quantile(values, q)
    end
    finite_extrema(A; default = (NaN, NaN)) = begin
        values = finite_values(A)
        isempty(values) ? default : extrema(values)
    end
    finite_sum(A) = sum(x -> isfinite(x) ? Float64(x) : 0.0, A)
    finite_maximum(A; default = NaN) = begin
        values = finite_values(A)
        isempty(values) ? default : maximum(values)
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
        valid = isfinite.(c.rho) .& isfinite.(c.vx) .& isfinite.(c.vy) .&
            isfinite.(c.vz) .& (c.rho .> 0)
        weights = ifelse.(valid, c.rho, 0.0)
        wsum = sum(weights)
        wsum > 0 || return (
            vbar = (NaN, NaN, NaN),
            dvx = fill(NaN, size(c.rho)), dvy = fill(NaN, size(c.rho)),
            dvz = fill(NaN, size(c.rho)),
            dv2 = fill(NaN, size(c.rho)),
        )
        vbar = (
            sum(weights .* ifelse.(valid, c.vx, 0.0)) / wsum,
            sum(weights .* ifelse.(valid, c.vy, 0.0)) / wsum,
            sum(weights .* ifelse.(valid, c.vz, 0.0)) / wsum,
        )
        dvx, dvy, dvz = c.vx .- vbar[1], c.vy .- vbar[2], c.vz .- vbar[3]
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

    Makie.get_tickvalues(::DecadeTicks, ::typeof(log10), vmin, vmax) =
        decade_tick_values(vmin, vmax)

    function Makie.get_ticks(::DecadeTicks, ::typeof(log10), formatter, vmin, vmax)
        values = decade_tick_values(vmin, vmax)
        values, latex_number.(values)
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
        axis = Axis(fig[1, 1],
            xlabel = L"N_{\mathrm H}\;[\mathrm{cm}^{-2}]",
            ylabel = ylabel, xscale = log10,
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

Each selected field creates one projected map for the active snapshot and line of sight. **Display** adds or removes a panel. **Logarithmic scale** transforms the displayed field before plotting; signed $B_{\mathrm{LOS}}$ values use a symmetric logarithm.

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

The energy maps show line-of-sight integrals, $\int E\,\mathrm{d}\ell$, and therefore represent energy per projected area. Arrows and line-integral convolution (LIC) trace the density-weighted plane-of-sky magnetic field. LIC follows the projected field periodically in both directions; its length, iterations, amplitude weighting, opacity, and deterministic seed are adjustable. **Contrast percentile** clips extreme values to preserve structure in the bulk of each map.
"""

# ╔═╡ 62440e86-b560-44ad-bb0a-43ae62e73fc3
md"""
---

## 12. Shared observational beam

The optional point-spread function is an elliptical Gaussian beam. Its major and minor FWHM can be specified in sky pixels or parsecs, and its position angle is measured counter-clockwise from the first displayed sky axis. The convolution is performed in Fourier space with periodic boundaries, consistently with the periodic simulation domain.

For polarized observations, the beam is applied to the linear Stokes quantities $I$, $Q$, and $U$ before computing

```math
P=\sqrt{Q^2+U^2},\qquad p=P/I.
```

This order captures beam depolarization and avoids the incorrect operation of smoothing $P$ or $p$ directly.

| Gaussian PSF setting | Control |
|:--|:--|
| Apply Gaussian beam | $(@bind apply_observational_beam PlutoUI.CheckBox(default = false)) |
| Beam-width unit | $(@bind observational_beam_unit PlutoUI.Select(["Sky pixels" => "pixel", "Parsecs" => "pc"]; default = "pixel")) |
| Major-axis FWHM | $(@bind observational_beam_fwhm_major PlutoUI.NumberField(0.1:0.1:1000.0; default = 3.0)) |
| Minor-axis FWHM | $(@bind observational_beam_fwhm_minor PlutoUI.NumberField(0.1:0.1:1000.0; default = 3.0)) |
| Position angle [$^\circ$] | $(@bind observational_beam_pa_deg PlutoUI.Slider(0.0:1.0:180.0; default = 0.0, show_value = true)) |
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

This section implements the line-of-sight geometry used by `GammaDust.jl`, `PhiDust.jl`, and `DustIQU.jl`. The magnetic inclination and polarization position angle satisfy $\gamma_{\mathrm d}=\arcsin(B_{\mathrm{LOS}}/|B|)$ and $\psi=\operatorname{atan2}(B_2,B_1)+\pi/2$. Optically thin emission follows $I_\nu,Q_\nu,U_\nu\propto B_\nu(T_{\mathrm d})\sigma_\nu n_{\mathrm H}\,\mathrm d\ell$, with the standard Stokes-$I$ correction $1-p_0(\cos^2\gamma_{\mathrm d}-2/3)$.

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
            ax = Axis(panel[1, 1], xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
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
            Colorbar(panel[1, 2], hm; label = as_latex(spec.label), tickformat = latex_ticklabels)
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
            axis = Axis(fig_dust_pixel_spectrum[1, index],
                xlabel = L"\nu\;[\mathrm{GHz}]", ylabel = spec.ylabel,
                xscale = log_dust_frequency_axis ? log10 : identity,
                xticks = log_dust_frequency_axis ? DECADE_TICKS : automatic,
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
    dust_stat_specs = Symbol[]
    show_dust_p_column && push!(dust_stat_specs, :column_relation)
    show_dust_p_pdf && push!(dust_stat_specs, :fraction_pdf)
    if isempty(dust_stat_specs)
        fig_dust_statistics = Figure(size = (900, 180))
        Label(fig_dust_statistics[1, 1], L"\mathrm{Select\ at\ least\ one\ dust\ statistic.}", fontsize = 20)
    else
        fig_dust_statistics = Figure(size = (520length(dust_stat_specs), 390))
        dust_N = vec(Float64.(column_density))
        dust_p = 100 .* vec(Float64.(dust_fraction))
        dust_valid = isfinite.(dust_N) .& isfinite.(dust_p) .& (dust_N .> 0) .& (dust_p .>= 0)
        dust_N_valid, dust_p_valid = dust_N[dust_valid], dust_p[dust_valid]
        for (index, statistic) in enumerate(dust_stat_specs)
            if statistic == :column_relation
                ax = Axis(fig_dust_statistics[1, index], xlabel = L"N_{\mathrm H}\;[\mathrm{cm}^{-2}]",
                    ylabel = L"100p_{\mathrm d}\;[\%]", xscale = log10,
                    xminorticks = IntervalsBetween(9), xminorticksvisible = true)
                sample_step = max(1, cld(length(dust_N_valid), 6000))
                sample = 1:sample_step:length(dust_N_valid)
                scatter!(ax, dust_N_valid[sample], dust_p_valid[sample];
                    color = (:dodgerblue3, 0.20), markersize = 4)
                if !isempty(dust_N_valid)
                    edges = 10 .^ range(extrema(log10.(dust_N_valid))...; length = 18)
                    centers, medians = Float64[], Float64[]
                    for bin in 1:length(edges)-1
                        members = (dust_N_valid .>= edges[bin]) .& (dust_N_valid .< edges[bin + 1])
                        any(members) || continue
                        push!(centers, sqrt(edges[bin] * edges[bin + 1]))
                        push!(medians, median(dust_p_valid[members]))
                    end
                    lines!(ax, centers, medians; color = MHD_COLORS[2], linewidth = 3)
                    scatter!(ax, centers, medians; color = MHD_COLORS[2], markersize = 7)
                end
            else
                ax = Axis(fig_dust_statistics[1, index], xlabel = L"100p_{\mathrm d}\;[\%]",
                    ylabel = L"\mathrm{PDF}(100p_{\mathrm d})")
                if !isempty(dust_p_valid)
                    upper = max(quantile(dust_p_valid, 0.999), 1.0e-6)
                    hist!(ax, clamp.(dust_p_valid, 0, upper); bins = range(0, upper; length = 45),
                        normalization = :pdf, color = (:seagreen3, 0.65), strokewidth = 1)
                end
            end
        end
    end
    fig_dust_p_column = polarization_column_figure(
        column_density, dust_fraction,
        L"100p_{\mathrm d}\;[\%]", MHD_COLORS[2])
    display_dust_statistics ? fig_dust_statistics : nothing
end

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
# ╟─a8ef96ab-0ddd-4eb2-a216-b7d96c2a9a08
# ╟─32110739-b60e-4592-856a-dd74f7a37401
# ╟─94a0a0dc-baf6-4e62-a51e-dc6124d98fd4
# ╟─c12d1f54-40b8-4865-9562-8dcb519f924a
# ╟─62440e86-b560-44ad-bb0a-43ae62e73fc3
# ╠═47b786d6-c7b5-44f4-946a-b8c485ad6380
# ╟─d6a2f4b1-59ac-4e77-a10a-4b74c0d89231
# ╟─130ccf03-d7cc-4b71-9210-dbc0e43dfa82
# ╟─5b3f6a91-246e-4cc3-8f68-164f7ff2f07c
# ╟─61c3b28c-9d62-4689-a85d-bc827e89641d
