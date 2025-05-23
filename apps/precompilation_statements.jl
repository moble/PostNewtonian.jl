using ZeroEccParamsFromPN: ZeroEccParamsFromPN

basic_args = ["--q=4.3", "--chiA=0.1,0.2,0.3", "--chiB=0.3,0.2,0.1"]

for arg ∈ ["--Omega0=0.015", "--D0=15"]
    args = [basic_args; arg]
    ZeroEccParamsFromPN.julia_main(args)
end
