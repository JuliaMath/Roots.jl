ENV["GKSwstype"] = "100"
using Documenter, DocumenterCodeBlocks
using DocumenterCitations, DocumenterLandingPage
using Roots


DocMeta.setdocmeta!(Roots, :DocTestSetup, :(using Roots); recursive=true)

makedocs(
    sitename = "Roots",
    format = Documenter.HTML(ansicolor=true),
    modules = [Roots],
    pages=[
        "Home" => "index.md",
        "Overview" => "roots.md",
        "Reference/API" => "reference.md",
        "Geometry" => "geometry-zero-finding.md"
    ],
    doctestfilters = [
        r"(?<=\d\.\d{12})\d+", # Ignore any digit after the 12th decimal place
    ],
plugins = [LandingPage(),
           CodeBlocks(),
           CitationBibliography("src/refs.bib")],
)

deploydocs(
    repo = "github.com/JuliaMath/Roots.jl"
)
