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
    ),
    (
        filename = "dust.jl",
        title = "Dust — thermal polarization",
        summary = "Synthetic thermal-dust Stokes emission, polarization maps, spectra, statistics, and polarization fraction versus column density.",
        start_heading = "## 13. Thermal-dust polarization",
        stop_heading = "## 14. Dichroic starlight polarization",
    ),
    (
        filename = "moose.jl",
        title = "MOOSE — Faraday tomography",
        summary = "Mock synchrotron emission, Faraday rotation, interferometric filtering, and RM synthesis.",
        start_heading = "## 16. MOOSE Faraday post-processing",
        stop_heading = "## 17. Polarization fraction versus intensity",
    ),
    (
        filename = "shine.jl",
        title = "SHINE — synthetic H I emission",
        summary = "21-cm radiative transfer, phase-separated H I columns, velocity moments, spectra, and RGB velocity composites.",
        start_heading = "## 18. SHINE H I post-processing",
        stop_heading = "## 19. Polarization fractions versus time",
    ),
    (
        filename = "zeeman.jl",
        title = "ZEEMAN — synthetic H I splitting",
        summary = "Synthetic Stokes I and V spectra, Zeeman field recovery, maps, and line-of-sight diagnostics.",
        start_heading = "## 15. H I Zeeman splitting",
        stop_heading = "## 16. MOOSE Faraday post-processing",
    ),
    (
        filename = "starlightpol.jl",
        title = "StarlightPol — dichroic polarization",
        summary = "Cell-by-cell Mueller propagation for background starlight through the simulated magnetized medium.",
        start_heading = "## 14. Dichroic starlight polarization",
        stop_heading = "## 15. H I Zeeman splitting",
    ),
]

const HERO_CELL_ID = Base.UUID("3ef88702-bef5-4eca-a151-df97aa7ec2c4")
const STYLE_CELL_ID = Base.UUID("8f1f1d9a-e6df-4dd1-b5b3-2d2c52c86686")
const LATEX_CELL_ID = Base.UUID("34d6b5f4-9c2d-42d5-9034-543aeb8ae151")
const TOC_CELL_ID = Base.UUID("bdc11245-d76c-45fe-b79d-7b64861f5f53")

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

function build_notebook(source, spec)
    cells = source.cells
    start_index = first_cell_containing(cells, spec.start_heading)
    stop_index = first_cell_containing(cells, spec.stop_heading)
    stop_index > start_index || error("Invalid section range for $(spec.filename)")

    seeds = cells[start_index:(stop_index - 1)]
    selected = dependency_closure(source.topology, seeds)

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
