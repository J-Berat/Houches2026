import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using Pluto

notebook_name = get(ENV, "DYNAMO_NOTEBOOK", "dynamo.jl")
notebook_path = abspath(joinpath(@__DIR__, notebook_name))
isfile(notebook_path) || error("Notebook not found: $notebook_path")
default_html_name = splitext(basename(notebook_path))[1] * ".html"
html_path = get(ENV, "DYNAMO_HTML_PATH", joinpath(@__DIR__, default_html_name))

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
            println(stderr, "\nPluto cell evaluation failed: ", cell.cell_id)
            show(stderr, MIME"text/plain"(), cell.output.body)
            println(stderr)
        end
        error("Cannot export notebook: $(length(failed_cells)) Pluto cell(s) failed.")
    end
    html = Pluto.generate_html(notebook; offline_bundle = true)
    write(html_path, html)
    println("HTML exported to: $html_path")
finally
    Pluto.SessionActions.shutdown(session, notebook; async = false)
end
