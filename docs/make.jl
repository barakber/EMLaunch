using Documenter

# Add src to load path to access EMLaunch module
push!(LOAD_PATH, "../src/")
using EMLaunch

makedocs(
    sitename = "EMLaunch.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [EMLaunch],
    pages = [
        "Home" => "index.md",
    ],
    remotes = nothing,  # Disable remote links for local builds
    checkdocs = :none   # Don't require all docstrings to be included
)

# Deploy documentation to gh-pages branch (only in CI)
if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/barakber/EMLaunch.git",
        devbranch = "main"
    )
end
