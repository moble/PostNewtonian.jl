using Pkg
Pkg.activate(@__DIR__)

using PackageCompiler
create_app(
    "ZeroEccParamsFromPN",
    "ZeroEccParamsFromPNApp";
    force=true,
    precompile_execution_file="precompilation_statements.jl",
)
