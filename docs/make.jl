ENV["GKSwstype"] = "100"
using Documenter
using Roots


DocMeta.setdocmeta!(Roots, :DocTestSetup, :(using Roots); recursive=true)

makedocs(
    sitename = "Roots",
    format = Documenter.HTML(ansicolor=true),
    modules = [Roots]
)

deploydocs(
    repo = "github.com/JuliaMath/Roots.jl"
)
