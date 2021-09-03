using Documenter
using Roots

DocMeta.setdocmeta!(Roots, :DocTestSetup, :(using Roots); recursive=true)

makedocs(
    sitename = "Roots",
    format = Documenter.HTML(),#ansicolor=true),
    modules = [Roots]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
