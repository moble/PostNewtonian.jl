"""
OLD DIRECTIONS:

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

@run_package_tests verbose = true

# ENV["JULIA_REFERENCETESTS_UPDATE"] = "true"
# @run_package_tests verbose = true filter=ti->(endswith(ti.filename, "reference_tests.jl") )
