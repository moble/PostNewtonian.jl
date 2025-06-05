using Pkg

# Ensure that the project inside the `ZeroEccParamsFromPN` directory is instantiated
# and working with the latest code.
Pkg.activate(joinpath(@__DIR__, "ZeroEccParamsFromPN"))
Pkg.develop(; path=dirname(@__DIR__))
Pkg.instantiate()

# Now switch to *this* project and ensure that it is instantiated.
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Create the app.  This may take ~5 minutes, and should create a bundle in this directory
# called `ZeroEccParamsFromPNApp`.
using PackageCompiler
create_app(
    "ZeroEccParamsFromPN",
    "ZeroEccParamsFromPNApp";
    force=true,  # Overwrite existing app
    precompile_execution_file="precompilation_statements.jl",
)
