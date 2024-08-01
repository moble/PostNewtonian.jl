using PostNewtonian
using Documenter

DocMeta.setdocmeta!(PostNewtonian, :DocTestSetup, :(using PostNewtonian); recursive=true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
    modules=[PostNewtonian],
    authors="Michael Boyle <michael.oliver.boyle@gmail.com> and contributors",
    repo="https://github.com/moble/PostNewtonian.jl/blob/{commit}{path}#{line}",
    sitename="PostNewtonian.jl",
    format=Documenter.HTML(; canonical="https://moble.github.io/PostNewtonian.jl"),
    pages=[
        "index.md"
        [
            file for file ∈ readdir(joinpath(@__DIR__, "src")) if
            file != "index.md" && splitext(file)[2] == ".md"
        ]
    ],
)

deploydocs(; repo="github.com/moble/PostNewtonian.jl")
