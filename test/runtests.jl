"""


A note about coverage: It's a well-hidden fact that you can turn coverage on and off by
adding certain comments around the code you don't want to measure:

    # COV_EXCL_START
    untested_code_that_wont_show_up_in_coverage()
    # COV_EXCL_STOP

"""

using TestItemRunner

function filter_files(filenames)
    # Check if the filename is in the list of filenames to allow
    return ti->any(endswith(ti.filename, f) for f ∈ filenames)
end

if length(ARGS) ≥ 1 && ARGS[1] == "update"
    ENV["JULIA_REFERENCETESTS_UPDATE"] = "true"
    @run_package_tests verbose = true filter=filter_files(["reference_tests.jl"])
elseif length(ARGS) > 0
    @run_package_tests verbose = true filter=filter_files(ARGS)
else
    @run_package_tests verbose = true
end
