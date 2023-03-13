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
        canonical="https://moble.github.io/PostNewtonian.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Interface" => [
            "High-level functions" => "interface.md",
            "GWFrames" => "gwframes.md",
            "Python" => "python.md",
            "Termination criteria" => "termination_criteria.md",
        ],
        "Internals" => [
            "Code structure" => "code_structure.md",
            "systems.md",
            "fundamental_variables.md",
            "derived_variables.md",
            "pn_expressions.md",
            "dynamics.md",
            "Evaluation" => "evaluation.md",
            "utilities.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/moble/PostNewtonian.jl",
    devbranch="main",
)
