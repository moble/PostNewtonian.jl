using PostNewtonian
using Documenter

DocMeta.setdocmeta!(
    PostNewtonian, :DocTestSetup, :(using PostNewtonian);
    recursive=true, warn=false
)

makedocs(;
    modules=[
        PostNewtonian,
        PostNewtonian.FundamentalVariables,
        PostNewtonian.DerivedVariables
    ],
    authors="Michael Boyle <michael.oliver.boyle@gmail.com> and contributors",
    repo = Remotes.GitHub("moble", "PostNewtonian.jl"),
    sitename="PostNewtonian.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://moble.github.io/PostNewtonian.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Interface" => [
            "High-level functions" => "interface/interface.md",
            "GWFrames" => "interface/gwframes.md",
            "Python" => "interface/python.md",
            "Symbolics" => "interface/symbolics.md",
            "Differentiation" => "interface/differentiation.md",
        ],
        "Internals" => [
            "Code structure" => "internals/code_structure.md",
            "internals/systems.md",
            "internals/fundamental_variables.md",
            "internals/derived_variables.md",
            "internals/pn_expressions.md",
            "internals/dynamics.md",
            "Waveforms" => "internals/waveforms.md",
            "internals/utilities.md",
        ],
        "Adding terms" => "adding_terms.md",
    ],
)

deploydocs(;
    repo="github.com/moble/PostNewtonian.jl",
    devbranch="main",
    push_preview=true,
)
