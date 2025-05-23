using Pkg
Pkg.activate(@__DIR__)

using PackageCompiler
create_app(
    "ZeroEccParamsFromPN",
    "ZeroEccParamsFromPNApp";
    force=true,  # Overwrite existing app
    precompile_execution_file="precompilation_statements.jl",
)
