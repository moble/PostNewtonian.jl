#!/usr/bin/env julia

# Run this script from any directory as
#
#   julia -t 4 scripts/docs.jl
#
# The docs will build and the browser should open automatically.  `LiveServer`
# will monitor the docs for any changes, then rebuild them and refresh the browser
# until this script is stopped.

using Dates: Dates
println("Building docs starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
#using Revise # Doesn't enable updates of docstrings in the output.
cd((@__DIR__) * "/..")
Pkg.activate("docs")
Pkg.instantiate()

using LiveServer
servedocs(; launch_browser=true)
