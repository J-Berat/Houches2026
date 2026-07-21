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
# Dynamo — MHD diagnostics

Spatial, statistical, temporal, and spectral diagnostics for the dynamo simulations.

> **Reactive mode.** Open `dynamo.jl` with Pluto. Selecting a repository, run, snapshot, or line of sight updates all dependent products automatically.

> **Lazy startup.** Run `run_pluto.jl` with `DYNAMO_NOTEBOOK=dynamo.jl`. Pluto starts without evaluating the expensive cells; run the result cells you need and Pluto will resolve their upstream dependencies.

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

# ╔═╡ 3df9ad9a-b865-41f5-8d7a-34a52d0292dd
Markdown.parse("""
### Repository status

| Item | Active value |
|:--|:--|
| Comparison parameter | $(comparison_parameter) |
| Resolved data root | **$(ROOT)** |
| Discovered runs | $(run_summary) |
""")

# ╔═╡ 290ecb0a-880a-44bd-9ddc-e6ee12b41a06
md"""
### Navigation and physical units

Choose the data cube and projection direction in the first table. Because the HDF5 files do not provide unit attributes, the second table defines the conversion from stored values to physical units. The default density scale is $10^{-12}\,\mathrm{g\,cm^{-3}}$ per stored unit; pressure, velocity, magnetic field, length, and time default to $\mathrm{erg\,cm^{-3}}$, $\mathrm{km\,s^{-1}}$, $\mathrm{G}$, $\mathrm{pc}$, and $\mathrm{Myr}$.

Comparative figures automatically detect the parameter varied from the selected path: Mach number for `VaryingMach`, grid resolution $N^3$ for `VaryingRes`, and $\chi=E_{\mathrm{comp}}/E_{\mathrm{sol}}$ for `VaryingRatio`.

#### FITS field convention

For a multi-extension FITS snapshot, use `EXTNAME` values such as `RHO`, `PRESSURE`, `VX`, `VY`, `VZ`, and either `BX`, `BY`, `BZ` or the face pairs `BX_L/BX_R`, `BY_L/BY_R`, `BZ_L/BZ_R`. The same names may be used as filenames inside one directory per snapshot. `L` and `TIME` are optional image extensions; the loader also accepts `LBOX`, `BOXSIZE`, `CDELT1..3`, `TIME`, `T`, or `SIMTIME` header keywords. If time metadata are absent, the final number in the filename or snapshot-directory name is used.

> **Analysis scope.** The **Run** and **Snapshot** controls select the active three-dimensional cube used by maps, PDFs, spectra, structure functions, and synthetic observations. **Simulations in comparative plots** independently selects the runs used by time histories, growth-rate relations, phase curves, and comparative energy or enstrophy diagnostics.
At startup, the active run is the only simulation selected for comparative plots. Add as many runs as needed with **Simulations in comparative plots**; every selected run is included.

- The adiabatic index $\gamma$ sets the sound-speed and thermal-energy convention.
- The mean particle mass $\mu m_{\mathrm H}$ defines $n=\rho/(\mu m_{\mathrm H})$ and $T=\mu m_{\mathrm H}P/(k_{\mathrm B}\rho)$.
- Magnetic energy and Alfvén speed use Gaussian CGS units: $E_{\mathrm B}=B^2/(8\pi)$ and $v_{\mathrm A}=B/\sqrt{4\pi\rho}$.
- **PDF weighting** selects equal-volume cell weighting or mass weighting.
- **Number of bins** controls the resolution of PDFs and density-binned diagnostics.
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
    active_time_value = run_times[selected_run][selected_snapshot] * time_unit_Myr
    active_time_text = isfinite(active_time_value) ?
        string(round(active_time_value; sigdigits = 6)) : "not available"
    comparison_runs_text = join(comparison_run_labels, ", ")
    open_cubes_text = join(analysis_series_labels, ", ")
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
    | Cubes opened | **$(length(analysis_series_labels))** |
    | Open cubes | **$(open_cubes_text)** |
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
            ax = Axis(panel[1, 1],
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
            Colorbar(panel[1, 2], hm, label = as_latex(spec.label), tickformat = latex_ticklabels)
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

The panels show number density $n$, magnetic-field strength $|B|$, and turbulent speed $|\delta\mathbf v|$ in physical units. Histograms are evaluated in $\log_{10}X$ bins. Their vertical axes are probability densities per dex and satisfy $\int (\mathrm{d}\mathcal P/\mathrm{d}\log_{10}X)\,\mathrm{d}\log_{10}X=1$.

**Display physical PDFs:** $(@bind display_pdfs PlutoUI.CheckBox(default = true))

| PDF panel | Display |
|:--|:--:|
| Number-density PDF | $(@bind show_pdf_density PlutoUI.CheckBox(default = true)) |
| Magnetic-field PDF | $(@bind show_pdf_magnetic PlutoUI.CheckBox(default = true)) |
| Turbulent-speed PDF | $(@bind show_pdf_velocity PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 496cbf2d-77a1-4a1a-b760-d4f8ea2ea9de
begin
    function density_pdf(values, weights, n)
        valid = isfinite.(values) .& isfinite.(weights) .& (weights .>= 0)
        values, weights = Float64.(values[valid]), Float64.(weights[valid])
        isempty(values) && return (Float64[], Float64[])
        lo, hi = quantile(values, (0.001, 0.999))
        if lo == hi
            lo, hi = lo - 0.5, hi + 0.5
        end
        edges = range(lo, hi; length = n + 1)
        h = fit(Histogram, values, Weights(weights), edges)
        centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
        pdf = h.weights ./ max(sum(h.weights .* diff(edges)), eps())
        centers, pdf
    end

    pdf_weights = pdf_weighting == "mass" ? vec(Float64.(cube.rho)) : ones(length(cube.rho))
    number_density_cells = number_density(cube.rho)
    magnetic_strength_uG = GAUSS_TO_MICROGAUSS .* mag.B
    turbulent_speed_kms = sqrt.(turb.dv2)
    mean_B = finite_mean(mag.B)
    sB = vec(safe_log10.(mag.B ./ mean_B))
    logn = vec(safe_log10.(number_density_cells))
    logBphysical = vec(safe_log10.(magnetic_strength_uG))
    logvphysical = vec(safe_log10.(turbulent_speed_kms))
    n_logx, n_p = density_pdf(logn, pdf_weights, nbins)
    B_logx, B_p = density_pdf(logBphysical, pdf_weights, nbins)
    v_logx, v_p = density_pdf(logvphysical, pdf_weights, nbins)
    logB_ratio_x, logB_ratio_p = density_pdf(sB, pdf_weights, nbins)
    n_x, B_x, v_x = 10.0 .^ n_logx, 10.0 .^ B_logx, 10.0 .^ v_logx
end

# ╔═╡ 1e0e6c0e-1ae0-40a1-a6f8-9c18fae91961
begin
    pdf_specs = NamedTuple[]
    show_pdf_density && push!(pdf_specs, (x = n_x, p = n_p,
        xlabel = L"n\;[\mathrm{cm}^{-3}]", color = MHD_COLORS[1]))
    show_pdf_magnetic && push!(pdf_specs, (x = B_x, p = B_p,
        xlabel = L"|B|\;[\mu\mathrm{G}]", color = MHD_COLORS[2]))
    show_pdf_velocity && push!(pdf_specs, (x = v_x, p = v_p,
        xlabel = L"|\delta v|\;[\mathrm{km\,s}^{-1}]", color = MHD_COLORS[3]))
    if isempty(pdf_specs)
        fig_pdf = Figure(size = (900, 180))
        Label(fig_pdf[1, 1], L"\mathrm{Select\ at\ least\ one\ PDF.}", fontsize = 20)
    else
        fig_pdf = Figure(size = (390length(pdf_specs), 360))
        for (j, spec) in enumerate(pdf_specs)
            ax = Axis(fig_pdf[1, j], xlabel = spec.xlabel,
                ylabel = L"\mathrm{d}\mathcal{P}/\mathrm{d}\log_{10}X", xscale = log10)
            stairs!(ax, spec.x, spec.p; color = spec.color, linewidth = 2.5, step = :center)
            ylims!(ax, low = 0)
        end
    end
    display_pdfs ? fig_pdf : nothing
end

# ╔═╡ fa155a62-da75-4530-b5e9-215fd4f66412
md"""
---

## 4. Thermodynamic phase diagram

This diagram shows the joint distribution of number density $n$ and thermal pressure $P/k_{\mathrm B}$ in the active three-dimensional cube. The plotted coordinates are $\log_{10}n$ and $\log_{10}(P/k_{\mathrm B})$; empty probability bins are masked.

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

    phase_logn = vec(safe_log10.(number_density(cube.rho)))
    phase_logpk = vec(safe_log10.(cube.P ./ K_B_CGS))
    phase_weights = pdf_weighting == "mass" ? vec(Float64.(cube.rho)) : ones(length(cube.rho))
    phase_data = phase_histogram(phase_logn, phase_logpk, phase_weights, phase_bins)
    phase_equilibrium = koyama_inutsuka_equilibrium()
end

# ╔═╡ 41b4eb12-889d-43b3-87c2-fc7cccf8679f
begin
    fig_phase = Figure(size = (680, 540))
    phase_axis = Axis(fig_phase[1, 1],
        xlabel = L"\log_{10}\!\left(n/\mathrm{cm}^{-3}\right)",
        ylabel = L"\log_{10}\!\left[(P/k_B)/(\mathrm{K\,cm}^{-3})\right]")
    phase_range = robust_colorrange(phase_data.log_probability, 99.0)
    phase_heatmap = heatmap!(phase_axis, phase_data.xcenters, phase_data.ycenters,
        phase_data.log_probability; colormap = :magma, colorrange = phase_range)
    if show_phase_equilibrium
        lines!(phase_axis, phase_equilibrium.logn, phase_equilibrium.logpk;
            color = :black, linewidth = 4.5)
        lines!(phase_axis, phase_equilibrium.logn, phase_equilibrium.logpk;
            color = :white, linewidth = 2.5, label = "n Λ(T) = Γ")
        axislegend(phase_axis; position = :rt, labelsize = 14)
    end
    phase_dx = length(phase_data.xcenters) > 1 ?
        phase_data.xcenters[2] - phase_data.xcenters[1] : 1.0
    phase_dy = length(phase_data.ycenters) > 1 ?
        phase_data.ycenters[2] - phase_data.ycenters[1] : 1.0
    xlims!(phase_axis, first(phase_data.xcenters) - phase_dx / 2,
        last(phase_data.xcenters) + phase_dx / 2)
    ylims!(phase_axis, first(phase_data.ycenters) - phase_dy / 2,
        last(phase_data.ycenters) + phase_dy / 2)
    Colorbar(fig_phase[1, 2], phase_heatmap,
        label = L"\log_{10}\mathcal{P}_{2\mathrm{D}}", tickformat = latex_ticklabels)
    colsize!(fig_phase.layout, 2, 22)
    display_phase_diagram ? fig_phase : nothing
end

# ╔═╡ 904ba663-d536-4b27-a379-4af927b0affb
begin
    function phase_magnetic_statistics(B, B2, rho, mask, weighting)
        valid = mask .& isfinite.(B) .& isfinite.(B2) .& isfinite.(rho)
        count(valid) == 0 && return (mean = NaN, rms = NaN)
        if weighting == "density"
            valid .&= rho .> 0
            weights = rho[valid]
            normalization = sum(weights)
            normalization > 0 || return (mean = NaN, rms = NaN)
            return (
                mean = GAUSS_TO_MICROGAUSS * sum(weights .* B[valid]) / normalization,
                rms = GAUSS_TO_MICROGAUSS * sqrt(sum(weights .* B2[valid]) / normalization),
            )
        end
        (
            mean = GAUSS_TO_MICROGAUSS * mean(B[valid]),
            rms = GAUSS_TO_MICROGAUSS * sqrt(mean(B2[valid])),
        )
    end

    function bulk_metrics(path, gamma)
        c = load_cube(path)
        m = magnetic_fields(c)
        u = turbulent_velocity(c)
        dynamic_valid = isfinite.(c.rho) .& isfinite.(u.dv2) .& (c.rho .> 0)
        dynamic_weights = ifelse.(dynamic_valid, c.rho, 0.0)
        mass = sum(dynamic_weights)
        vrms = mass > 0 ? sqrt(sum(dynamic_weights .* ifelse.(dynamic_valid, u.dv2, 0.0)) / mass) : NaN
        cs2_kms2 = gamma .* c.P ./ c.rho ./ KM_CM^2
        va2_kms2 = m.B2 ./ (4pi .* c.rho) ./ KM_CM^2
        cs_valid = dynamic_valid .& isfinite.(cs2_kms2) .& (cs2_kms2 .>= 0)
        va_valid = dynamic_valid .& isfinite.(va2_kms2) .& (va2_kms2 .>= 0)
        cs_mass = sum(ifelse.(cs_valid, c.rho, 0.0))
        va_mass = sum(ifelse.(va_valid, c.rho, 0.0))
        cs_rms = cs_mass > 0 ? sqrt(sum(ifelse.(cs_valid, c.rho .* cs2_kms2, 0.0)) / cs_mass) : NaN
        va_rms = va_mass > 0 ? sqrt(sum(ifelse.(va_valid, c.rho .* va2_kms2, 0.0)) / va_mass) : NaN
        Ekin = finite_sum(0.5 .* c.rho .* u.dv2 .* KM_CM^2)
        Emag = finite_sum(m.B2 ./ (8pi))
        Etherm = finite_sum(gamma > 1 + sqrt(eps(Float64)) ? c.P ./ (gamma - 1) : c.P)
        (
            t = c.t,
            mach = vrms / max(cs_rms, eps()),
            mach_alfven = vrms / max(va_rms, eps()),
            Bmean = GAUSS_TO_MICROGAUSS * finite_mean(m.B),
            Brms = GAUSS_TO_MICROGAUSS * sqrt(finite_mean(m.B2)),
            energy_ratio = Ekin > 0 ? Emag / Ekin : NaN,
            kin_mag = Ekin > 0 && Emag > 0 ? Ekin / Emag : NaN,
            therm_mag = Etherm > 0 && Emag > 0 ? Etherm / Emag : NaN,
            kin_therm = Ekin > 0 && Etherm > 0 ? Ekin / Etherm : NaN,
        )
    end

    function phase_magnetic_metrics(path, molecular_weight, cold_boundary_K,
            warm_boundary_K, phase_weighting)
        c = load_cube(path)
        m = magnetic_fields(c)
        temperature = molecular_weight * M_H_CGS .* c.P ./ (K_B_CGS .* c.rho)
        cold_boundary = min(cold_boundary_K, warm_boundary_K)
        warm_boundary = max(cold_boundary_K, warm_boundary_K)
        cold_B = phase_magnetic_statistics(m.B, m.B2, c.rho,
            temperature .< cold_boundary, phase_weighting)
        lukewarm_B = phase_magnetic_statistics(m.B, m.B2, c.rho,
            (temperature .>= cold_boundary) .& (temperature .< warm_boundary), phase_weighting)
        warm_B = phase_magnetic_statistics(m.B, m.B2, c.rho,
            temperature .>= warm_boundary, phase_weighting)
        (
            t = c.t,
            B_cold_mean = cold_B.mean,
            B_cold_rms = cold_B.rms,
            B_lukewarm_mean = lukewarm_B.mean,
            B_lukewarm_rms = lukewarm_B.rms,
            B_warm_mean = warm_B.mean,
            B_warm_rms = warm_B.rms,
        )
    end

    all_series = Dict(
        label => bulk_metrics.(run_files[label], gamma)
        for label in analysis_series_labels
    )

    phase_B_series_by_run = display_phase_B_time ? Dict(
        label => phase_magnetic_metrics.(run_files[label],
            Float64(mean_molecular_weight), Float64(phase_cold_boundary_K),
            Float64(phase_warm_boundary_K), phase_B_weighting)
        for label in comparison_run_labels
    ) : Dict{String, Any}()
end

# ╔═╡ 36aef377-3de7-435a-af83-3a90421e3159
md"""
---

## 5. Magnetic amplification and growth-rate fit

The fitted physical model is $B(t)=A\exp[\Gamma_B(t-t_0)]$, where $t_0$ is the first valid snapshot time and both $A$ and $\Gamma_B$ are inferred from the snapshots in the **Fit window**. Equivalently, the regression is $\ln B=a+\Gamma_B(t-t_0)$ with a free intercept. The shaded interval marks the selected fit range. This avoids treating the first magnetic-field measurement as exact or giving it special leverage in the regression.

The figures show $\ln(B/B_0)$ on linear axes, with $B_0$ the first valid measured field used only as a plotting normalization. The reported $R^2$ uses the usual mean-centred total sum of squares appropriate to this intercept fit. The theoretical comparison curves remain anchored at the first valid measured snapshot. Because $E_{\mathrm B}\propto B^2$, magnetic energy grows at rate $2\Gamma_B$.

**Display the magnetic growth-rate fit:** $(@bind display_growth_fit PlutoUI.CheckBox(default = true))

| Setting | Control |
|:--|:--|
| Fitted field | $(@bind growth_fit_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Fit window (snapshot indices) | $(@bind growth_fit_window PlutoUI.RangeSlider(1:length(all_series[selected_run]); default = min(2, length(all_series[selected_run])):min(4, length(all_series[selected_run])), show_value = true)) |
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

# ╔═╡ 71c8ea26-d2ad-4430-9265-0b28d45bba1c
md"""
### Interval magnetic-growth rate

The local magnetic-growth rate is evaluated between consecutive snapshots rather than with one fit over the complete selected window:

```math
\Gamma_{B,i+1/2}
=\frac{\ln B_{i+1}-\ln B_i}{t_{i+1}-t_i}.
```

Time, sonic Mach number, and Alfvénic Mach number are assigned their midpoint values over the same interval. Positive $\Gamma_B$ indicates magnetic amplification; negative $\Gamma_B$ indicates decay.

**Display interval growth-rate relations:** $(@bind display_gamma_relations PlutoUI.CheckBox(default = true))

| Growth-rate setting | Control |
|:--|:--|
| Magnetic field | $(@bind gamma_relation_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Simulations | Global multi-selection in section 1 |
| $\Gamma_B(t)$ | $(@bind show_gamma_time PlutoUI.CheckBox(default = true)) |
| $\Gamma_B(\mathcal{M})$ | $(@bind show_gamma_mach PlutoUI.CheckBox(default = true)) |
| $\Gamma_B(\mathcal{M}_{\mathrm A})$ | $(@bind show_gamma_alfven_mach PlutoUI.CheckBox(default = true)) |
"""

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
            gamma_axis = Axis(fig_gamma_relations[1, panel_index],
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

# ╔═╡ cd7e037c-62f5-4682-92c2-92af7169d692
md"""
### Multi-snapshot normalized magnetic evolution

Selected snapshots are superimposed in the same figure using

```math
\ln\!\left(\frac{B_i}{B_0}\right),
```

where $B_0$ is the first available snapshot of each run. Color identifies the run and marker shape identifies the snapshot. All axes remain linear because the logarithm is applied to the magnetic-field ratio itself.

**Display multi-snapshot normalized magnetic evolution:** $(@bind display_normalized_B_relations PlutoUI.CheckBox(default = true))

| Normalized-field relation | Control |
|:--|:--|
| Magnetic field | $(@bind normalized_B_relation_field PlutoUI.Select(["Mean field ⟨B⟩", "RMS field Bᵣₘₛ"]; default = "Mean field ⟨B⟩")) |
| Snapshots displayed | $(@bind normalized_B_snapshot_indices PlutoUI.MultiSelect(collect(1:maximum_snapshot_count); default = unique([1, cld(maximum_snapshot_count, 2), maximum_snapshot_count]))) |
| Simulations displayed | Global multi-selection in section 1 |
| $\ln(B/B_0)$ versus $t$ | $(@bind show_normalized_B_time PlutoUI.CheckBox(default = true)) |
| $\ln(B/B_0)$ versus $\mathcal{M}$ | $(@bind show_normalized_B_mach PlutoUI.CheckBox(default = true)) |
| $\ln(B/B_0)$ versus $\mathcal{M}_{\mathrm A}$ | $(@bind show_normalized_B_alfven_mach PlutoUI.CheckBox(default = true)) |
"""

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
    normalized_B_marker_cycle = [:circle, :rect, :utriangle, :diamond, :cross,
        :xcross, :dtriangle, :pentagon, :hexagon]
    normalized_B_markers = Dict(
        snapshot => normalized_B_marker_cycle[mod1(position, length(normalized_B_marker_cycle))]
        for (position, snapshot) in enumerate(normalized_B_valid_snapshots)
    )

    if isempty(normalized_B_panel_specs) || isempty(normalized_B_valid_snapshots) ||
            isempty(normalized_B_runs)
        fig_normalized_B_relations = Figure(size = (900, 180))
        Label(fig_normalized_B_relations[1, 1],
            L"\mathrm{Select\ at\ least\ one\ relation\ and\ one\ snapshot.}", fontsize = 20)
    else
        normalized_B_ncols = length(normalized_B_panel_specs)
        fig_normalized_B_relations = Figure(size = (430normalized_B_ncols, 500))
        for (panel_index, spec) in enumerate(normalized_B_panel_specs)
            axis = Axis(fig_normalized_B_relations[1, panel_index],
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
                for (x, y, snapshot) in zip(horizontal[valid], logarithmic_ratio[valid], available[valid])
                    scatter!(axis, [x], [y]; color = run_colors[label],
                        marker = normalized_B_markers[snapshot], markersize = 11)
                end
            end
        end
        run_elements = [LineElement(color = run_colors[label], linewidth = 2.4)
            for label in normalized_B_runs]
        snapshot_elements = [MarkerElement(marker = normalized_B_markers[snapshot],
            color = :gray25, markersize = 11) for snapshot in normalized_B_valid_snapshots]
        Legend(fig_normalized_B_relations[2, 1:length(normalized_B_panel_specs)],
            [run_elements, snapshot_elements],
            [legend_run_label.(normalized_B_runs),
                ["i = $(snapshot)" for snapshot in normalized_B_valid_snapshots]],
            ["Run", "Snapshot"];
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
        fig_growth = Figure(size = (550growth_panel_count, 430))
        growth_column = 0
        if show_growth_fit_panel
            growth_column += 1
            ag1 = Axis(fig_growth[1, growth_column], xlabel = L"t\;[\mathrm{Myr}]",
                ylabel = L"\ln(B/B_0)")
            lines!(ag1, growth_times, ln_normalized_B; color = run_colors[selected_run],
                linewidth = 2.5, label = "Data")
            scatter!(ag1, growth_times, ln_normalized_B; color = run_colors[selected_run], markersize = 8)
            if show_growth_fit && isfinite(growth_fit.gamma)
                fit_t = growth_times[growth_indices]
                fit_log_ratio = growth_fit.log_amplitude - log(growth_B0) .+
                    growth_fit.gamma .* growth_elapsed_time[growth_indices]
                lines!(ag1, fit_t, fit_log_ratio; color = :black, linewidth = 3,
                    linestyle = :dash,
                    label = string("ln(A/B₀) + ΓB(t − t₀), ",
                        legend_rate_label(growth_fit.gamma; fitted = true)))
            end
            growth_has_interval && vspan!(ag1,
                growth_times[growth_first_index], growth_times[growth_last_index];
                color = (:gray50, 0.10))
            axislegend(ag1, position = :lt)
        end
        if show_growth_theory_panel
            growth_column += 1
            ag2 = Axis(fig_growth[1, growth_column], xlabel = L"t\;[\mathrm{Myr}]",
                ylabel = L"\ln(B/B_0)")
            lines!(ag2, growth_times, ln_normalized_B; color = run_colors[selected_run],
                linewidth = 2.5, label = "Data")
            scatter!(ag2, growth_times, ln_normalized_B; color = run_colors[selected_run], markersize = 7)
            for (Γ, curve, color) in zip(theory_gammas, theory_B_curves, theory_colors)
                lines!(ag2, growth_times, log.(curve ./ growth_B0); color, linewidth = 2,
                    linestyle = :dot,
                    label = legend_rate_label(Γ))
            end
            if show_growth_fit && isfinite(growth_fit.gamma)
                lines!(ag2, growth_times, log.(fitted_ratio_curve); color = :black,
                    linewidth = 2.5, linestyle = :dashdot,
                    label = "ln(A/B₀) + ΓB, fit (t − t₀)")
            end
            growth_has_interval && vspan!(ag2,
                growth_times[growth_first_index], growth_times[growth_last_index];
                color = (:gray50, 0.10))
            axislegend(ag2, position = :lt)
        end
    end
    display_growth_fit ? fig_growth : nothing
end

# ╔═╡ dcc8f8f3-daaf-4d4b-92a3-4919ed5e36de
md"""
---

## 6. Global time evolution

Selected panels compare every discovered run as a function of physical time. Turbulent velocity is measured after subtracting the mass-weighted bulk motion.

$\mathcal{M}=v_{\mathrm{rms}}/c_{s,\mathrm{rms}}$ measures compressibility, while $\mathcal{M}_{\mathrm A}=v_{\mathrm{rms}}/v_{\mathrm A,\mathrm{rms}}$ compares turbulence with Alfvén-wave propagation. Values above unity are supersonic or super-Alfvénic, respectively. Dotted magnetic-field curves show the theoretical exponentials selected in section 5.

**Display global time evolution:** $(@bind display_global_evolution PlutoUI.CheckBox(default = true))

| Time-evolution panel | Display |
|:--|:--:|
| Sonic Mach number | $(@bind show_time_mach PlutoUI.CheckBox(default = true)) |
| Alfvénic Mach number | $(@bind show_time_alfven PlutoUI.CheckBox(default = true)) |
| Magnetic field | $(@bind show_time_magnetic PlutoUI.CheckBox(default = true)) |
| Magnetic-to-kinetic energy ratio | $(@bind show_time_energy_ratio PlutoUI.CheckBox(default = true)) |
"""

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
                Axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]", ylabel = L"\mathcal{M}") :
                key == :alfven ?
                Axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]", ylabel = L"\mathcal{M}_{\mathrm{A}}") :
                key == :magnetic ?
                Axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]",
                    ylabel = L"B\;[\mu\mathrm{G}]", yscale = log_B_time ? log10 : identity) :
                Axis(fig_time[row, col], xlabel = L"t\;[\mathrm{Myr}]",
                    ylabel = L"E_B/E_{\mathrm{kin}}", yscale = log10)
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

**Display magnetic field by phase:** $(@bind display_phase_B_time PlutoUI.CheckBox(default = true))

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

# ╔═╡ d56a8ba3-aa42-4351-a065-87275032a342
begin
    phase_B_specs = NamedTuple[]
    show_phase_B_cold && push!(phase_B_specs,
        (key = :cold, label = "CNM", linestyle = :solid, marker = :circle))
    show_phase_B_lukewarm && push!(phase_B_specs,
        (key = :lukewarm, label = "LNM", linestyle = :dash, marker = :rect))
    show_phase_B_warm && push!(phase_B_specs,
        (key = :warm, label = "WNM", linestyle = :dot, marker = :utriangle))

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
        phase_B_axis = Axis(fig_phase_B_time[1, 1],
            xlabel = L"t\;[\mathrm{Myr}]", ylabel = phase_B_ylabel,
            yscale = log_phase_B_time ? log10 : identity)
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
                color = run_colors[label], marker = spec.marker, markersize = 7)
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

The map shows the line-of-sight mean of the transformed field $\log_{10}(B/\langle B\rangle)$. The PDF uses every cell in the active cube. A value of zero corresponds to the cube-mean field strength; positive and negative values identify locally stronger and weaker fields.

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
            axmap = Axis(map_panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            limB = max(maximum(abs, filter(isfinite, vec(logB_map)); init = 0.0),
                sqrt(eps(Float64)))
            hb = heatmap!(axmap, sky_coordinates[1], sky_coordinates[2], logB_map;
                colormap = :balance, colorrange = (-limB, limB))
            Colorbar(map_panel[1, 2], hb,
                label = L"\log_{10}(B/\langle B\rangle)", tickformat = latex_ticklabels)
            colsize!(map_panel, 2, 22)
        end
        if show_logB_histogram
            logB_column += 1
            axhist = Axis(fig_logB[1, logB_column],
                xlabel = L"\log_{10}(B/\langle B\rangle)", ylabel = L"\mathcal{P}")
            lines!(axhist, logB_ratio_x, logB_ratio_p;
                color = MHD_COLORS[2], linewidth = 2.5)
        end
    end
    display_normalized_field ? fig_logB : nothing
end

# ╔═╡ 298bd579-bb28-48b7-8c55-ea74804b9837
md"""
---

## 8. Energy partition by density and time

### Density-binned energy ratios

The figure compares three energy ratios as functions of physical number density. Color identifies the simulation run and line style identifies the snapshot. The horizontal reference at unity marks equal energies; values above unity indicate that the numerator dominates.

Each point is a ratio of energies summed within one density bin, rather than a mean of cell-by-cell ratios. The CGS energy densities are $E_{\mathrm{kin}}=\rho|\delta\mathbf v|^2/2$, $E_{\mathrm{mag}}=B^2/(8\pi)$, and $E_{\mathrm{therm}}=P/(\gamma-1)$ in $\mathrm{erg\,cm}^{-3}$. For $\gamma=1$, the notebook adopts the isothermal convention $E_{\mathrm{therm}}=P$. Cubes selected for the comparative profiles are read and reduced to these per-cell quantities in a separate Pluto cell; changing $\gamma$ or the density bins therefore rebins the cached arrays without rereading the RAMSES files.

**Display energy ratios by density:** $(@bind display_energy_ratios PlutoUI.CheckBox(default = true))

**Display energy ratios versus time:** $(@bind display_energy_time PlutoUI.CheckBox(default = true))

| Setting | Control |
|:--|:--|
| Snapshots shown in energy reports | $(@bind energy_snapshot_indices PlutoUI.MultiSelect(collect(1:maximum_snapshot_count); default = unique([min(2, maximum_snapshot_count), cld(maximum_snapshot_count, 2), maximum_snapshot_count]))) |
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
    energy_profile_samples = Dict(
        (label, snapshot_index) => energy_density_samples(
            load_cube(run_files[label][snapshot_index]))
        for label in comparison_run_labels for snapshot_index in valid_energy_snapshots
        if snapshot_index <= length(run_files[label])
    )
end

# ╔═╡ 2c4762f5-4780-4bfa-a320-e4102d9f5d6c
begin
    energy_profiles = Dict(
        (label, snapshot_index) => energy_ratios_by_density(
            energy_profile_samples[(label, snapshot_index)], density_edges, gamma,
        )
        for (label, snapshot_index) in keys(energy_profile_samples)
    )
end

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
            energy_axes[key] = Axis(fig_energy[1, index],
                xlabel = L"n\;[\mathrm{cm}^{-3}]", ylabel = ylabel,
                xscale = log10, yscale = log10)
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
            ax = Axis(fig_energy_time[1, index],
                xlabel = L"t\;[\mathrm{Myr}]", ylabel = spec.ylabel, yscale = log10)
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
| Density–family-parameter heatmap | $(@bind show_enstrophy_parameter_heatmap PlutoUI.CheckBox(default = true)) |
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
            axis = Axis(panel[1, 1],
                xlabel = latexstring(sky_labels[1], "/\\mathrm{pc}"),
                ylabel = latexstring(sky_labels[2], "/\\mathrm{pc}"))
            heat = heatmap!(axis, sky_coordinates[1], sky_coordinates[2], spec.data;
                colormap = spec.colormap)
            Colorbar(panel[1, 2], heat; label = as_latex(spec.label), tickformat = latex_ticklabels)
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
    enstrophy_profiles = Dict(
        label => enstrophy_by_density(
            load_cube(run_files[label][enstrophy_snapshot_indices[label]]),
            density_edges, enstrophy_density_weighting,
        ) for label in comparison_run_labels
    )
    enstrophy_profile_matrix = hcat([
        enstrophy_profiles[label] for label in comparison_run_labels]...)

    function enstrophy_parameter_label(label)
        if comparison_kind == :mach
            index = enstrophy_snapshot_indices[label]
            mach_value = all_series[label][index].mach
            return string("𝓜 = ", round(mach_value; sigdigits = 3))
        end
        legend_run_label(label)
    end

    function enstrophy_parameter_ticklabel(label)
        if comparison_kind == :mach
            index = enstrophy_snapshot_indices[label]
            mach_value = all_series[label][index].mach
            return latexstring("\\mathcal{M}=", round(mach_value; sigdigits = 3))
        end
        latex_run_label(label)
    end

    enstrophy_parameter_axis_label = comparison_kind == :resolution ?
        L"N^3" : comparison_kind == :ratio ?
        L"\chi=E_{\mathrm{comp}}/E_{\mathrm{sol}}" : comparison_kind == :mach ?
        L"\mathcal{M}" : L"\mathrm{simulation}"
    # Keep categorical tick labels as atomic strings. LaTeXString labels can be
    # split into a different number of Makie text blocks while Pluto updates the
    # corresponding tick positions, triggering a ComputePipeline length error.
    enstrophy_parameter_ticklabels = String.(
        enstrophy_parameter_ticklabel.(comparison_run_labels))
end

# ╔═╡ 873f7ef2-719b-4ae6-b015-1a23c6c27836
begin
    enstrophy_density_panels = Symbol[]
    show_enstrophy_density_profiles && push!(enstrophy_density_panels, :profiles)
    show_enstrophy_parameter_heatmap && push!(enstrophy_density_panels, :parameter_heatmap)
    if isempty(enstrophy_density_panels)
        fig_enstrophy_density = Figure(size = (900, 180))
        Label(fig_enstrophy_density[1, 1],
            L"\mathrm{Select\ at\ least\ one\ density-binned\ enstrophy\ diagnostic.}",
            fontsize = 20)
    else
        fig_enstrophy_density = Figure(size = (610length(enstrophy_density_panels), 500))
        for (panel_index, panel_kind) in enumerate(enstrophy_density_panels)
            if panel_kind == :profiles
                axis = Axis(fig_enstrophy_density[1, panel_index],
                    xlabel = L"n\;[\mathrm{cm}^{-3}]",
                    ylabel = L"\langle\mathcal{E}_{\omega}\rangle_n\;[\mathrm{Myr}^{-2}]",
                    xscale = log10, yscale = log10,
                    xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                    xminorticksvisible = true, yminorticksvisible = true)
                for label in comparison_run_labels
                    profile = enstrophy_profiles[label]
                    valid = isfinite.(profile) .& (profile .> 0)
                    lines!(axis, density_number_centers[valid], profile[valid];
                        color = run_colors[label], linewidth = 2.5, label = enstrophy_parameter_label(label))
                    scatter!(axis, density_number_centers[valid], profile[valid];
                        color = run_colors[label], markersize = 5)
                end
                axislegend(axis; position = :lt, framevisible = false)
            else
                panel = fig_enstrophy_density[1, panel_index] = GridLayout()
                parameter_positions = collect(1:length(comparison_run_labels))
                axis = Axis(panel[1, 1], xlabel = enstrophy_parameter_axis_label,
                    ylabel = L"n\;[\mathrm{cm}^{-3}]", yscale = log10,
                    yminorticks = IntervalsBetween(9), yminorticksvisible = true,
                    xticks = parameter_positions,
                    xtickformat = values -> [
                        enstrophy_parameter_ticklabels[clamp(round(Int, value), 1,
                            length(enstrophy_parameter_ticklabels))] for value in values
                    ])
                log_enstrophy = safe_log10.(enstrophy_profile_matrix')
                heat = heatmap!(axis, parameter_positions, density_number_centers,
                    log_enstrophy; colormap = :magma)
                Colorbar(panel[1, 2], heat;
                    label = L"\log_{10}\!\left(\langle\mathcal{E}_{\omega}\rangle_n/\mathrm{Myr}^{-2}\right)",
                    tickformat = latex_ticklabels)
                colsize!(panel, 2, 22)
            end
        end
    end
    display_enstrophy_density ? fig_enstrophy_density : nothing
end

# ╔═╡ a8558c31-7dcf-433e-9950-a59e9acf158b
md"""
---

## 10. Isotropic power spectra

The panels show number-density, turbulent-velocity, and vorticity spectra on base-10 logarithmic axes in both $k$ and spectral power. Fourier power is summed in spherical shells and normalized consistently with Parseval's theorem. Dividing shell power by $\Delta k=2\pi/L$ gives spectral densities satisfying $\int E_f(k)\,\mathrm{d}k=\langle|f|^2\rangle$, with the stated velocity prefactor. Velocity and vorticity include all three vector components, and $k$ is expressed in $\mathrm{pc}^{-1}$.

The first shells contain few Fourier modes and are correspondingly noisy. The shaded region $k<3\,\Delta k$ marks these box-scale modes; it should be excluded when estimating an inertial-range slope. The threshold is a visual reliability guide, not a claim that an inertial range necessarily begins at $3\,\Delta k$.

**Display density, velocity, and vorticity power spectra:** $(@bind display_power_spectra PlutoUI.CheckBox(default = true))

| Power-spectrum panel | Display |
|:--|:--:|
| Number-density spectrum | $(@bind show_spectrum_density PlutoUI.CheckBox(default = true)) |
| Velocity spectrum | $(@bind show_spectrum_velocity PlutoUI.CheckBox(default = true)) |
| Vorticity spectrum | $(@bind show_spectrum_vorticity PlutoUI.CheckBox(default = true)) |
"""

# ╔═╡ 3a731972-3404-478c-a572-00a05ab652b1
begin
    spectrum_specs = NamedTuple[]
    show_spectrum_density && push!(spectrum_specs,
        (k = krho, power = Prho, ylabel = L"E_n(k)\;[\mathrm{cm}^{-6}\,\mathrm{pc}]",
            color = MHD_COLORS[1]))
    show_spectrum_velocity && push!(spectrum_specs,
        (k = kv, power = Pv,
            ylabel = L"E_v(k)\;[(\mathrm{km\,s}^{-1})^2\,\mathrm{pc}]",
            color = MHD_COLORS[3]))
    show_spectrum_vorticity && push!(spectrum_specs,
        (k = komega, power = Pomega,
            ylabel = L"E_\omega(k)\;[\mathrm{Myr}^{-2}\,\mathrm{pc}]",
            color = MHD_COLORS[5]))
    if isempty(spectrum_specs)
        fig_spectra = Figure(size = (900, 180))
        Label(fig_spectra[1, 1], L"\mathrm{Select\ at\ least\ one\ power\ spectrum.}", fontsize = 20)
    else
        fig_spectra = Figure(size = (400length(spectrum_specs), 390))
        for (index, spec) in enumerate(spectrum_specs)
            axis = Axis(fig_spectra[1, index], xlabel = L"k\;[\mathrm{pc}^{-1}]",
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
            lines!(axis, spec.k[valid], spec.power[valid]; color = spec.color, linewidth = 2.5)
            scatter!(axis, spec.k[valid], spec.power[valid]; color = spec.color, markersize = 5)
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
            axis = Axis(fig_structure[1, index], xlabel = L"\ell\;[\mathrm{pc}]",
                ylabel = spec.ylabel, xscale = log10, yscale = log10,
                xticks = DECADE_TICKS, yticks = DECADE_TICKS,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                xminorticksvisible = true, yminorticksvisible = true)
            lines!(axis, structure_separations_pc, spec.values;
                color = spec.color, linewidth = 2.5)
            scatter!(axis, structure_separations_pc, spec.values;
                color = spec.color, markersize = 6)
        end
    end
    display_structure_functions ? fig_structure : nothing
end

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
