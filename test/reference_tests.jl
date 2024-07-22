@testitem "Reference tests" begin

    using Quaternionic
    using PostNewtonian
    using ReferenceTests
    using HDF5
    using FileIO
    using Random

    function approx(reference, actual; atol=1e-8, rtol=1e-8)
        for (k, v1) in reference
            if haskey(actual, k)
                v2 = actual[k]
                if size(v1)≠size(v2) || !isapprox(v1, v2; atol=atol, rtol=rtol)
                    @warn "Mismatch for key $k:\n\t$v1\n\t≉\n\t$v2"
                    return false
                end
            else
                @warn """Key "$k" found in "reference" file but not in "actual"."""
                return false
            end
        end
        for (k, v2) in actual
            if !haskey(reference, k)
                @warn """Key "$k" found in "actual" file but not in "reference"."""
                return false
            end
        end
        return true
    end

    function compare_v̇(file_name, args...; kwargs...)
        file_name = joinpath("reference_values", file_name)
        rng = Random.Xoshiro(1234);
        pnsystems = [rand(rng, BBH) for _ ∈ 1:1_000]
        @test_reference file_name Dict(
            "2PNOrder" => (@. Int(2pn_order(pnsystems))),
            "state" => reduce(hcat, getproperty.(pnsystems, :state)),
            "v̇_numerator" => PostNewtonian.v̇_numerator.(pnsystems),
            "v̇_denominator" => PostNewtonian.v̇_denominator.(pnsystems),
        ) by=(r,a)->approx(r,a,atol=1e-12,rtol=1e-12)
    end

    function compare_inspiral(file_name, args...; kwargs...)
        file_name = joinpath("reference_values", file_name)
        inspiral = orbital_evolution(args...; kwargs...)
        pnsystem = try
            inspiral.prob.p
        catch
            inspiral.p
        end
        @test_reference file_name Dict(
            "2PNOrder" => Int(2pn_order(pnsystem)),
            "state" => pnsystem.state,
            "t" => inspiral.t,
            "u" => reduce(hcat, inspiral.u)
        ) by=approx
    end

    function compare_waveforms(file_name, inspiral, f; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=typemax(Int))
        file_name = joinpath("reference_values", file_name)
        h = f(inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
        @test_reference file_name Dict(
            "h" => h
        ) by=approx
    end

    M₁ = 0.8
    M₂ = 0.2
    χ⃗₁ = QuatVec(0.1, 0.2, 0.3)
    χ⃗₂ = QuatVec(0.6, 0.4, 0.1)
    Ωᵢ = PostNewtonian.Ω(v=0.3, M=M₁+M₂)
    Ωₑ = PostNewtonian.Ω(v=0.4, M=M₁+M₂)

    inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Ωₑ)

    compare_v̇("vdot.h5")

    # These tests rely too heavily on ODE integration being consistent, which seems to be
    # highly architecture-dependent, so we disable them on CI.
    @show get(ENV, "CI", "not present!")
    if get(ENV, "CI", "false") == "false"
        compare_inspiral("inspiral1.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ)
        compare_inspiral("inspiral2.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Ω₁=Ωᵢ/2)
        compare_inspiral("inspiral3.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, saveat=3.7)
        compare_inspiral("inspiral4.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, saveat=inspiral.t[3:7:end-9])
        compare_inspiral("inspiral5.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Ωₑ)
        compare_inspiral("inspiral6.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, saves_per_orbit=3)
        compare_inspiral("inspiral7.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, PNOrder=1//2)
        compare_inspiral("inspiral8.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, PNOrder=3//2)
        compare_inspiral("inspiral9.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, PNOrder=5//2)
        compare_inspiral("inspiral10.h5", M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, PNOrder=7//2)

        compare_waveforms("waveform1.h5", inspiral, coorbital_waveform)
        compare_waveforms("waveform2.h5", inspiral, coorbital_waveform, ℓₘᵢₙ=1, ℓₘₐₓ=4)
        compare_waveforms("waveform3.h5", inspiral, coorbital_waveform, PNOrder=1//2)
        compare_waveforms("waveform4.h5", inspiral, coorbital_waveform, PNOrder=3//2)
        compare_waveforms("waveform5.h5", inspiral, coorbital_waveform, PNOrder=5//2)
        compare_waveforms("waveform6.h5", inspiral, coorbital_waveform, PNOrder=7//2)
        compare_waveforms("waveform7.h5", inspiral, inertial_waveform)
        compare_waveforms("waveform8.h5", inspiral, inertial_waveform, ℓₘᵢₙ=1, ℓₘₐₓ=4)
        compare_waveforms("waveform9.h5", inspiral, inertial_waveform, PNOrder=1//2)
        compare_waveforms("waveform10.h5", inspiral, inertial_waveform, PNOrder=3//2)
        compare_waveforms("waveform11.h5", inspiral, inertial_waveform, PNOrder=5//2)
        compare_waveforms("waveform12.h5", inspiral, inertial_waveform, PNOrder=7//2)
    end

end
