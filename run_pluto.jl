import Pkg

Pkg.activate(@__DIR__)

if VERSION < v"1.11" || VERSION >= v"1.13"
    @warn "Ce dépôt est testé avec Julia 1.11 et 1.12." julia_version = VERSION
end

try
    @eval using Pluto
catch
    println(stderr, "\nL'environnement Julia n'est pas encore installé.")
    println(stderr, "Exécuter une fois :")
    println(stderr, "  julia --project=. -e 'import Pkg; Pkg.instantiate()'\n")
    rethrow()
end

if !isdefined(Pluto, :DEFAULT_PRECEDENCE_HEURISTIC)
    @eval Pluto const DEFAULT_PRECEDENCE_HEURISTIC =
        PlutoDependencyExplorer.DEFAULT_PRECEDENCE_HEURISTIC
end

const AVAILABLE_NOTEBOOKS = [
    ("Complet — toutes les analyses", "dynamo_diagnostics.jl"),
    ("Dynamo — diagnostics MHD", joinpath("notebooks", "dynamo.jl")),
    ("Dust — polarisation thermique", joinpath("notebooks", "dust.jl")),
    ("MOOSE — tomographie Faraday", joinpath("notebooks", "moose.jl")),
    ("SHINE — émission H I", joinpath("notebooks", "shine.jl")),
    ("ZEEMAN — dédoublement H I", joinpath("notebooks", "zeeman.jl")),
    (
        "StarlightPol — polarisation dichroïque",
        joinpath("notebooks", "starlightpol.jl"),
    ),
]

function resolve_notebook(requested)
    cleaned = lowercase(strip(requested))
    isempty(cleaned) && return nothing

    numeric = tryparse(Int, cleaned)
    if !isnothing(numeric) && numeric in eachindex(AVAILABLE_NOTEBOOKS)
        return AVAILABLE_NOTEBOOKS[numeric][2]
    end

    requested_name =
        endswith(cleaned, ".jl") ? basename(cleaned) : basename(cleaned) * ".jl"
    match_index = findfirst(AVAILABLE_NOTEBOOKS) do entry
        lowercase(basename(entry[2])) == requested_name
    end
    isnothing(match_index) ? nothing : AVAILABLE_NOTEBOOKS[match_index][2]
end

function ask_for_notebook()
    println("\nNotebooks disponibles :\n")
    for (index, (label, path)) in enumerate(AVAILABLE_NOTEBOOKS)
        println("  ", index, ". ", label, "  [", path, "]")
    end

    while true
        print("\nNotebook à ouvrir [1-", length(AVAILABLE_NOTEBOOKS), "] : ")
        flush(stdout)
        eof(stdin) && error("Aucun notebook sélectionné.")
        notebook = resolve_notebook(readline())
        !isnothing(notebook) && return notebook
        println("Choix invalide.")
    end
end

requested = isempty(ARGS) ? get(ENV, "DYNAMO_NOTEBOOK", "") : first(ARGS)
relative_notebook = isempty(strip(requested)) ?
    ask_for_notebook() : resolve_notebook(requested)
isnothing(relative_notebook) && error("Notebook inconnu : $requested")

notebook = abspath(joinpath(@__DIR__, relative_notebook))
isfile(notebook) || error(
    "Notebook introuvable : $notebook\n" *
    "Régénérer les notebooks avec : " *
    "julia --project=. tools/split_notebooks.jl",
)

parse_boolean(value) =
    lowercase(strip(value)) in ("1", "true", "yes", "y", "on")

pluto_host = get(ENV, "PLUTO_HOST", "127.0.0.1")
pluto_port = parse(Int, get(ENV, "PLUTO_PORT", "1234"))
launch_browser = parse_boolean(get(
    ENV,
    "PLUTO_LAUNCH_BROWSER",
    string(Sys.isapple() || Sys.iswindows()),
))

println("Notebook : ", notebook)
println("Démarrage paresseux : aucun cube n'est ouvert automatiquement.")
Pluto.run(
    notebook_path_suggestion = notebook,
    host = pluto_host,
    port = pluto_port,
    launch_browser = launch_browser,
    run_notebook_on_load = false,
)
