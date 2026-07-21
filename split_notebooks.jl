import Pkg

Pkg.activate(@__DIR__)

using Pluto

const SOURCE_NOTEBOOK = joinpath(@__DIR__, "dynamo_diagnostics.jl")

const NOTEBOOK_SPECS = [
    (
        filename = "dynamo.jl",
        title = "Dynamo — MHD diagnostics",
        summary = "Spatial, statistical, temporal, and spectral diagnostics for the dynamo simulations.",
        start_heading = "# Interactive MHD diagnostics",
        stop_heading = "## 13. Thermal-dust polarization",
        figures = [
            ("heatmaps", "Projected heatmaps", :fig_maps),
            ("pdfs", "Probability density functions", :fig_pdf),
            ("phase_diagram", "Pressure-density phase diagram", :fig_phase),
            ("time_evolution", "Global time evolution", :fig_time),
            ("phase_magnetic_time", "Magnetic field by thermal phase", :fig_phase_B_time),
            ("magnetic_fit", "Magnetic exponential fit", :fig_growth),
            ("growth_rate_relations", "Growth rate versus time and Mach numbers", :fig_gamma_relations),
            ("normalized_magnetic_relations", "Normalized magnetic evolution", :fig_normalized_B_relations),
            ("normalized_magnetic_field", "Normalized magnetic-field distribution", :fig_logB),
            ("energy_ratios", "Energy ratios by density", :fig_energy),
            ("energy_time", "Energy ratios versus time", :fig_energy_time),
            ("vorticity", "Vorticity map", :fig_vorticity),
            ("enstrophy_density", "Enstrophy by density", :fig_enstrophy_density),
            ("power_spectra", "Power spectra", :fig_spectra),
            ("structure_functions", "Structure functions", :fig_structure),
        ],
    ),
    (
        filename = "dust.jl",
        title = "Dust — thermal polarization",
        summary = "Synthetic thermal-dust Stokes emission, polarization maps, spectra, statistics, and polarization fraction versus column density.",
        start_heading = "## 13. Thermal-dust polarization",
        stop_heading = "## 14. Dichroic starlight polarization",
        figures = [
            ("dust_polarization", "Dust-polarization maps", :fig_dust),
            ("dust_structure", "Dust observable structure functions", :fig_dust_structure),
            ("dust_pixel_spectrum", "Dust Stokes spectrum", :fig_dust_pixel_spectrum),
            ("dust_statistics", "Dust comparative statistics", :fig_dust_statistics),
            ("dust_p_column", "Dust polarization versus column density", :fig_dust_p_column),
        ],
    ),
    (
        filename = "moose.jl",
        title = "MOOSE — Faraday tomography",
        summary = "Mock synchrotron emission, Faraday rotation, interferometric filtering, and RM synthesis.",
        start_heading = "## 16. MOOSE Faraday post-processing",
        stop_heading = "## 17. Polarization fraction versus intensity",
        figures = [
            ("moose", "MOOSE maps", :fig_moose),
            ("moose_structure", "MOOSE observable structure functions", :fig_moose_structure),
            ("moose_tomography", "MOOSE Faraday tomography", :fig_moose_tomography),
            ("moose_p_column", "Faraday polarization versus column density", :fig_moose_p_column),
        ],
    ),
    (
        filename = "shine.jl",
        title = "SHINE — synthetic H I emission",
        summary = "21-cm radiative transfer, phase-separated H I columns, velocity moments, spectra, and RGB velocity composites.",
        start_heading = "## 18. SHINE H I post-processing",
        stop_heading = "## 19. Polarization fractions versus time",
        figures = [
            ("shine", "SHINE H I maps", :fig_shine),
            ("shine_structure", "SHINE observable structure functions", :fig_shine_structure),
            ("shine_rgb", "SHINE velocity RGB composite", :fig_shine_rgb),
            ("shine_spectrum", "SHINE H I spectrum", :fig_shine_spectrum),
        ],
    ),
    (
        filename = "zeeman.jl",
        title = "ZEEMAN — synthetic H I splitting",
        summary = "Synthetic Stokes I and V spectra, Zeeman field recovery, maps, and line-of-sight diagnostics.",
        start_heading = "## 15. H I Zeeman splitting",
        stop_heading = "## 16. MOOSE Faraday post-processing",
        figures = [
            ("zeeman_maps", "Zeeman maps", :fig_zeeman_maps),
            ("zeeman_structure", "Zeeman observable structure functions", :fig_zeeman_structure),
            ("zeeman_spectra", "Zeeman Stokes spectra", :fig_zeeman_spectra),
            ("zeeman_p_column", "Zeeman polarization versus column density", :fig_zeeman_p_column),
        ],
    ),
    (
        filename = "starlightpol.jl",
        title = "StarlightPol — dichroic polarization",
        summary = "Cell-by-cell Mueller propagation for background starlight through the simulated magnetized medium.",
        start_heading = "## 14. Dichroic starlight polarization",
        stop_heading = "## 15. H I Zeeman splitting",
        figures = [
            ("starlight_maps", "Starlight-polarization maps", :fig_starlight_maps),
            ("starlight_structure", "Starlight observable structure functions", :fig_starlight_structure),
            ("starlight_profiles", "Starlight sight-line profiles", :fig_starlight_profiles),
            ("starlight_p_column", "Starlight polarization versus column density", :fig_starlight_p_column),
        ],
    ),
]

const HERO_CELL_ID = Base.UUID("3ef88702-bef5-4eca-a151-df97aa7ec2c4")
const STYLE_CELL_ID = Base.UUID("8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686")
const LATEX_CELL_ID = Base.UUID("34d6b5f4-9c2d-42d5-9034-543aeb8ae151")
const TOC_CELL_ID = Base.UUID("bdc11245-d76c-45fe-b79d-7b64861f5f53")
const EXPORT_OPTIONS_CELL_ID = Base.UUID("e1000001-6f8c-4d0c-9a10-000000000001")
const EXPORT_CONTROLS_CELL_ID = Base.UUID("e1000002-6f8c-4d0c-9a10-000000000002")
const EXPORT_DOWNLOAD_CELL_ID = Base.UUID("e1000003-6f8c-4d0c-9a10-000000000003")

function first_cell_containing(cells, text)
    index = findfirst(cell -> occursin(text, cell.code), cells)
    isnothing(index) && error("Heading not found: $text")
    index
end

function dependency_closure(topology, seed_cells)
    selected = Set(seed_cells)
    pending = collect(seed_cells)
    while !isempty(pending)
        cell = pop!(pending)
        references = topology.nodes[cell].references
        for upstream in Pluto.where_assigned(topology, references)
            if upstream ∉ selected
                push!(selected, upstream)
                push!(pending, upstream)
            end
        end
    end
    selected
end

function hero_code(title, summary, filename)
    """md\"\"\"
    # $title

    $summary

    > **Reactive mode.** Open `$filename` with Pluto. Selecting a repository, run, snapshot, or line of sight updates all dependent products automatically.

    > **Lazy startup.** Run `run_pluto.jl` with `DYNAMO_NOTEBOOK=$filename`. Pluto starts without evaluating the expensive cells; run the result cells you need and Pluto will resolve their upstream dependencies.

    All dimensional quantities are converted to the physical units shown on their axes or colorbars. Projected means are density weighted unless stated otherwise, and periodic boundaries are used for spatial operations.
    \"\"\""""
end

function export_options_code(figures)
    options = join(["        \"$key\" => \"$label\"," for (key, label, _) in figures], "\n")
    registry = join(["        \"$key\" => $figure," for (key, _, figure) in figures], "\n")
    """begin
    export_figure_options = [
$options
    ]
    export_figure_registry = Dict(
$registry
    )
    nothing
end"""
end

function export_controls_code(figures)
    default_key = first(figures)[1]
    string(
        "md\"\"\"\n---\n\n## Figure export\n\n",
        "| Export setting | Control |\n|:--|:--|\n",
        "| Figure | \$(@bind export_figure_key PlutoUI.Select(export_figure_options; default = \"$default_key\")) |\n",
        "| Format | \$(@bind export_figure_format PlutoUI.Select([\"PNG\", \"PDF\"]; default = \"PDF\")) |\n",
        "\"\"\"",
    )
end


function export_download_code(filename)
    notebook_slug = splitext(filename)[1]
    """begin
    export_extension = lowercase(export_figure_format)
    export_mime = export_figure_format == \"PNG\" ? MIME\"image/png\"() : MIME\"application/pdf\"()
    export_buffer = IOBuffer()
    show(export_buffer, export_mime, export_figure_registry[export_figure_key])
    export_bytes = take!(export_buffer)
    export_run_slug = replace(lowercase(selected_run), r\"[^a-z0-9]+\" => \"_\")
    export_filename = \"$(notebook_slug)_\$(export_figure_key)_\$(export_run_slug)_snapshot_\$(lpad(selected_snapshot, 3, '0')).\$(export_extension)\"
    PlutoUI.DownloadButton(export_bytes, export_filename)
end"""
end

function build_notebook(source, spec)
    cells = source.cells
    start_index = first_cell_containing(cells, spec.start_heading)
    stop_index = first_cell_containing(cells, spec.stop_heading)
    stop_index > start_index || error("Invalid section range for $(spec.filename)")

    seeds = cells[start_index:(stop_index - 1)]
    selected = dependency_closure(source.topology, seeds)

    if spec.filename == "dynamo.jl"
        beam_start = first_cell_containing(cells, "## 12. Shared observational beam")
        for cell in cells[beam_start:(stop_index - 1)]
            delete!(selected, cell)
        end
    end

    for required_id in (STYLE_CELL_ID, LATEX_CELL_ID, TOC_CELL_ID)
        push!(selected, source.cells_dict[required_id])
    end

    output_cells = Pluto.Cell[]
    for cell in cells
        if cell.cell_id == HERO_CELL_ID
            push!(output_cells, Pluto.Cell(HERO_CELL_ID,
                hero_code(spec.title, spec.summary, spec.filename)))
        elseif cell in selected
            push!(output_cells, cell)
        end
    end

    any(cell -> cell.cell_id == HERO_CELL_ID, output_cells) ||
        insert!(output_cells, 3, Pluto.Cell(HERO_CELL_ID,
            hero_code(spec.title, spec.summary, spec.filename)))

    push!(output_cells,
        Pluto.Cell(EXPORT_OPTIONS_CELL_ID, export_options_code(spec.figures)),
        Pluto.Cell(EXPORT_CONTROLS_CELL_ID, export_controls_code(spec.figures)),
        Pluto.Cell(EXPORT_DOWNLOAD_CELL_ID, export_download_code(spec.filename)),
    )

    output_path = joinpath(@__DIR__, spec.filename)
    output_notebook = Pluto.Notebook(output_cells, output_path)
    output_notebook.nbpkg_ctx = source.nbpkg_ctx
    output_notebook.nbpkg_ctx_instantiated = source.nbpkg_ctx_instantiated
    Pluto.save_notebook(output_notebook, output_path)
    println(rpad(spec.filename, 20), length(output_cells), " cells")
end

source = Pluto.load_notebook_nobackup(SOURCE_NOTEBOOK)
source.topology = Pluto.updated_topology(source.topology, source, source.cells)
foreach(spec -> build_notebook(source, spec), NOTEBOOK_SPECS)
