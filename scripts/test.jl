#!/usr/bin/env julia

# Run this script (from any directory) either as an executable with `test.jl`, or with
# arguments to julia as in `julia --threads=auto test.jl`.  In either case, you can follow
# the command with any number of file names to restrict to running tests only in those
# files.  Alternatively, you can add the argument `update` to update the reference tests
# (and skip all others).

using Dates: Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd(dirname(@__DIR__))
Pkg.activate(".")
Pkg.instantiate()

try
    Δt = @elapsed Pkg.test("PostNewtonian"; coverage=true, test_args=ARGS)
    println("Running tests took $Δt seconds.")
catch e
    println("Tests failed; proceeding to coverage")
end

Pkg.activate()  # Activate Julia's base (home) directory
using Coverage
cd(dirname(@__DIR__))
coverage = Coverage.process_folder()
Coverage.writefile("lcov.info", coverage)
Coverage.clean_folder(".")
