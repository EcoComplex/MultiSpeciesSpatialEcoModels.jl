using Documenter
using MultiSpeciesSpatialEcoModels

makedocs(
    sitename = "MultiSpeciesSpatialEcoModels.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/EcoComplex/MultiSpeciesSpatialEcoModels.jl.git",
    devbranch = "main"
)