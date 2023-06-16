# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, qColoringCardinality

makedocs(
    modules = [qColoringCardinality],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Moe Kayali",
    sitename = "qColoringCardinality.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/mkyl/qColoringCardinality.jl.git",
    push_preview = true
)
