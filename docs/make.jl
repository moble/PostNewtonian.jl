using PostNewtonian
using Documenter

DocMeta.setdocmeta!(PostNewtonian, :DocTestSetup, :(using PostNewtonian); recursive=true)

makedocs(;
    modules=[PostNewtonian],
    authors="Michael Boyle <michael.oliver.boyle@gmail.com> and contributors",
    repo="https://github.com/moble/PostNewtonian.jl/blob/{commit}{path}#{line}",
    sitename="PostNewtonian.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://moble.github.io/PostNewtonian.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/moble/PostNewtonian.jl",
    devbranch="main",
)