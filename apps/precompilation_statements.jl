using ZeroEccParamsFromPN: ZeroEccParamsFromPN

# We want to precompile all branches of the code, so we have to input a variety of
# arguments.  We can look at the `ArgParse` stuff in the main function to see what options
# are available.
#
# The first group of arguments describes the mass and spin parameters of the binary system.
# These shouldn't affect the branches significantly, so we can just use one set of those
# arguments.
#
# The next group requires exactly one of `D0`, `Omega0, `tMerger`, or `NOrbits` to be
# passed.  Again, the precise values shouldn't affect branching much, so we use just one set
# of arguments.
#
# The next group is not yet implemented, so we skip it for now.
#
# The final group includes two flags.  The first allows us to skip the up-down-instability
# check, which would only skip a branch, so we just leave it out.  The second is for the D0
# hack, so we run both with and without to exercise both options.

basic_args = ["--q=4.3", "--chiA=0.1,0.2,0.3", "--chiB=0.3,0.2,0.1"]

for arg ∈ ["--Omega0=0.015", "--D0=15", "--tMerger=4000", "--NOrbits=20"]
    args = [basic_args; arg]
    ZeroEccParamsFromPN.julia_main(args)
    ZeroEccParamsFromPN.julia_main([args; "--experimental_D0_hack"])
end
