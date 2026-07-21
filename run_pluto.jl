import Pkg

Pkg.activate(@__DIR__)

if VERSION < v"1.11" || VERSION >= v"1.13"
    @warn "This repository is tested with Julia 1.11 and 1.12." julia_version = VERSION
end

active_manifest = Pkg.Types.Context().env.manifest_file
println("Julia version: ", VERSION)
println("Dependency manifest: ", basename(active_manifest))

try
    @eval using Pluto
catch
    println(stderr, "\nThe Julia environment is not installed yet.")
    println(stderr, "Run this command once, then start the launcher again:")
    println(stderr, "  julia --project=. -e 'import Pkg; Pkg.instantiate()'\n")
    rethrow()
end

# Pluto 1.0.3 references this dependency constant without its module prefix
# when lazy startup is enabled. Define the missing alias only for affected
# versions; later Pluto releases that define it themselves are left unchanged.
if !isdefined(Pluto, :DEFAULT_PRECEDENCE_HEURISTIC)
    @eval Pluto const DEFAULT_PRECEDENCE_HEURISTIC =
        PlutoDependencyExplorer.DEFAULT_PRECEDENCE_HEURISTIC
end

const AVAILABLE_NOTEBOOKS = [
    ("Dynamo — MHD diagnostics", "dynamo.jl"),
    ("Dust — thermal polarization", "dust.jl"),
    ("MOOSE — Faraday tomography", "moose.jl"),
    ("SHINE — synthetic H I emission", "shine.jl"),
    ("ZEEMAN — synthetic H I splitting", "zeeman.jl"),
    ("StarlightPol — dichroic polarization", "starlightpol.jl"),
]

function resolve_notebook_choice(choice)
    cleaned = lowercase(strip(choice))
    isempty(cleaned) && return nothing

    numeric_choice = tryparse(Int, cleaned)
    if !isnothing(numeric_choice) && numeric_choice in eachindex(AVAILABLE_NOTEBOOKS)
        return AVAILABLE_NOTEBOOKS[numeric_choice][2]
    end

    normalized_name = endswith(cleaned, ".jl") ? cleaned : cleaned * ".jl"
    match_index = findfirst(entry -> lowercase(entry[2]) == normalized_name,
        AVAILABLE_NOTEBOOKS)
    isnothing(match_index) ? nothing : AVAILABLE_NOTEBOOKS[match_index][2]
end

function ask_for_notebook()
    println("\nAvailable Pluto notebooks:\n")
    for (index, (label, filename)) in enumerate(AVAILABLE_NOTEBOOKS)
        println("  $index. $label  [$filename]")
    end

    while true
        print("\nWhich notebook do you want to open? [1-$(length(AVAILABLE_NOTEBOOKS))]: ")
        flush(stdout)
        eof(stdin) && error("No notebook choice was provided.")
        notebook_name = resolve_notebook_choice(readline())
        !isnothing(notebook_name) && return notebook_name
        println("Invalid choice. Enter a number from 1 to $(length(AVAILABLE_NOTEBOOKS)), or a notebook name.")
    end
end

requested_notebook = isempty(ARGS) ? get(ENV, "DYNAMO_NOTEBOOK", "") : first(ARGS)
notebook_name = if isempty(strip(requested_notebook))
    ask_for_notebook()
else
    resolved = resolve_notebook_choice(requested_notebook)
    isnothing(resolved) && error("Unknown notebook: $requested_notebook")
    resolved
end

notebook = abspath(joinpath(@__DIR__, notebook_name))
isfile(notebook) || error("Notebook not found: $notebook")

parse_boolean(value) = lowercase(strip(value)) in ("1", "true", "yes", "y", "on")
pluto_host = get(ENV, "HOST", "127.0.0.1")
pluto_port = parse(Int, get(ENV, "PORT", "1234"))
default_browser_setting = Sys.isapple() || Sys.iswindows()
launch_browser = parse_boolean(get(ENV, "LAUNCH_BROWSER",
    string(default_browser_setting)))

println("Opening Pluto for: $notebook")
println("Lazy startup enabled: select a cell and run it with Shift+Enter.")
!launch_browser && println("Open http://$pluto_host:$pluto_port in your browser (use an SSH tunnel for a remote server).")
Pluto.run(
    notebook_path_suggestion = notebook,
    host = pluto_host,
    port = pluto_port,
    launch_browser = launch_browser,
    run_notebook_on_load = false,
)
