"""
Run this script (from this directory) as

    time julia --code-coverage=tracefile-%p.info --code-coverage=user --project=. ./runtests.jl

Then, if you have `lcov` installed, you should also have `genhtml`, and you can run this

    genhtml tracefile-<your_PID>.info --output-directory coverage/ && open coverage/index.html

to view the coverage locally as HTML.  I find that this sometimes requires
removing files that aren't really there from the .info file.

It's a well-hidden fact that you can turn coverage on and off by adding certain comments around the
code you don't want to measure:

    # COV_EXCL_START
    untested_code_that_wont_show_up_in_coverage()
    # COV_EXCL_STOP

"""

using TestItemRunner

@run_package_tests verbose=true

# using PostNewtonian
# using Test

# using Logging
# using Random
# using SciMLBase
# using OrdinaryDiffEq
# using Symbolics
# using SymbolicUtils
# using DoubleFloats
# using Quaternionic

# import Aqua
# import ExplicitImports

# enabled_tests = lowercase.(ARGS)

# help = ("help" ∈ enabled_tests || "--help" ∈ enabled_tests)
# helptests = []

# # This block is cribbed from StaticArrays.jl/test/runtests.jl
# function addtests(fname)
#     key = lowercase(splitext(fname)[1])
#     if help
#         push!(helptests, key)
#     else
#         if isempty(enabled_tests) || key in enabled_tests
#             println("Running $key.jl")
#             Random.seed!(42)
#             include(fname)
#         end
#     end
# end

# @testset verbose=true "PostNewtonian" begin
#     addtests("aqua.jl")
#     addtests("macros.jl")
#     addtests("binding_energy.jl")
#     addtests("up_down_instability.jl")
#     addtests("orbital_evolution.jl")
#     addtests("gwframes.jl")
#     addtests("truncated_series_monoid.jl")
#     addtests("pnexpansion.jl")
# end

# if help
#     println()
#     println("Pass no args to run all tests, or select one or more of the following:")
#     for helptest in helptests
#         println("    ", helptest)
#     end
# end
