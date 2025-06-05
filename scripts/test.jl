#!/usr/bin/env julia

using Dates: Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
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
cd((@__DIR__) * "/..")
coverage = Coverage.process_folder()
Coverage.writefile("lcov.info", coverage)
Coverage.clean_folder(".")
