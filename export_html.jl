import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using Pluto

requested_notebook =
    get(ENV, "DYNAMO_NOTEBOOK", "dynamo_diagnostics.jl")
candidate_paths = [
    abspath(joinpath(@__DIR__, requested_notebook)),
    abspath(joinpath(@__DIR__, "notebooks", basename(requested_notebook))),
]
notebook_index = findfirst(isfile, candidate_paths)
isnothing(notebook_index) &&
    error("Notebook introuvable : $requested_notebook")
notebook_path = candidate_paths[notebook_index]
html_path = abspath(get(
    ENV,
    "DYNAMO_HTML_PATH",
    joinpath(
        @__DIR__,
        "exports",
        splitext(basename(notebook_path))[1] * ".html",
    ),
))
mkpath(dirname(html_path))

session = Pluto.ServerSession()
session.options.evaluation.workspace_use_distributed = false
session.options.evaluation.run_notebook_on_load = true
session.options.server.disable_writing_notebook_files = true

notebook = Pluto.SessionActions.open(
    session,
    notebook_path;
    run_async = false,
)

try
    failed_cells = filter(cell -> cell.errored, notebook.cells)
    if !isempty(failed_cells)
        for cell in failed_cells
            println(stderr, "\nÉchec de la cellule Pluto : ", cell.cell_id)
            show(stderr, MIME"text/plain"(), cell.output.body)
            println(stderr)
        end
        error("Export interrompu : $(length(failed_cells)) cellule(s) en erreur.")
    end
    html = Pluto.generate_html(notebook; offline_bundle = true)
    write(html_path, html)
    println("HTML enregistré : ", html_path)
finally
    Pluto.SessionActions.shutdown(session, notebook; async = false)
end
