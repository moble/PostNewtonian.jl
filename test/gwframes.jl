@testset verbose=true "GWFrames" begin
    Random.seed!(1234)
    T = Float64
    M₁ = T(5//8)
    M₂ = T(3//8)
    χ⃗₁ = normalize(randn(QuatVec{T})) * rand(T(0):T(1//1_000_000):T(1))
    χ⃗₂ = normalize(randn(QuatVec{T})) * rand(T(0):T(1//1_000_000):T(1))
    Ωᵢ = T(1//64)
    Ω₁ = Ωᵢ/2
    Ωₑ = 1
    Rᵢ = exp(normalize(randn(QuatVec{T})) * rand(T(0):T(1//1_000_000):T(1//1_000)))

    Mₜₒₜ = M₁+M₂
    q = M₁/M₂
    χ⃗ₛ = χₛ(M₁, M₂, χ⃗₁, χ⃗₂)
    χ⃗ₐ = χₐ(M₁, M₂, χ⃗₁, χ⃗₂)
    vᵢ = v(Ω=Ωᵢ,M=M₁+M₂)
    v₁ = v(Ω=Ω₁,M=M₁+M₂)
    vₑ = min(v(Ω=Ωₑ, M=M₁+M₂), 1)

    uᵢ = [M₁; M₂; χ⃗₁.vec; χ⃗₂.vec; Rᵢ.components; vᵢ]

    Approximant = "TaylorT1"
    delta = δ(M₁, M₂)
    chi1_i = χ⃗₁.vec
    chi2_i = χ⃗₂.vec
    Omega_orb_i = Ωᵢ
    Omega_orb_0 = Ω₁

    # Just make sure we don't get any errors
    GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i)
    GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, Omega_orb_0)
    GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, Omega_orb_0, [Rᵢ.components...])
    GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, Omega_orb_0, [Rᵢ.components...], 37)

    # Test to ensure we don't get info with default value of `quiet`, and we *do* when `quiet=false`
    @test_logs min_level=Logging.Info GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, Omega_orb_0)
    @test_logs min_level=Logging.Info GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0)
    @test_logs min_level=Logging.Info GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0, quiet=true)
    @test_logs (:info,r"Terminating forwards") (:info,r"Terminating backwards") GWFrames.PNWaveform(
        Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0, quiet=false)

    @testset verbose=true "MinStepsPerOrbit" begin
        for MinStepsPerOrbit ∈ [17, 31, 32, 33, 64]
            w = GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0, MinStepsPerOrbit)
            for Φᵢ ∈ w.Phi[begin:end-MinStepsPerOrbit]
                NSteps = sum(@. Φᵢ ≤ w.Phi ≤ Φᵢ + 2π * (1 - 1/100MinStepsPerOrbit))
                @test NSteps == MinStepsPerOrbit
            end
        end
    end

    @testset verbose=true "Uniform time steps" begin
        for dt in [1, 0.3, 0.1]
            w = GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0, dt)
            difft = diff(w.t)
            @test dt-minimum(difft) < length(w.t)*10eps(dt)
            @test maximum(difft)-dt < length(w.t)*10eps(dt)
        end
    end

end
