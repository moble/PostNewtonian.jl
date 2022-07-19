@testset verbose=true "Inspiral" begin
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

    forwards_termination = (
        "Terminating forwards evolution because the PN parameter 𝑣 "
        * "has reached 𝑣ₑ=$(vₑ).  This is ideal."
    )
    backwards_termination = (
        "Terminating backwards evolution because the PN parameter 𝑣 "
        * "has reached 𝑣₁=$(v₁).  This is ideal."
    )

    # Check for termination info
    sol1 = @test_logs (:info,forwards_termination) inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Rᵢ=Rᵢ)
    sol2 = @test_logs (:info,forwards_termination) (:info,backwards_termination) inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Ω₁=Ωᵢ/2, Rᵢ=Rᵢ)

    # Check endpoint values
    @test sol1.retcode == :Terminated
    @test sol1[end, 1] == vᵢ
    @test sol1[1] ≈ uᵢ
    @test sol1[end, end] ≈ vₑ

    @test sol2.retcode == :Terminated
    @test sol2[end, 1] ≈ v₁
    iᵢ = argmin(abs.(sol2.t .- 0.0))  # Assuming uᵢ corresponds to t==0.0
    @test sol2[iᵢ] ≈ uᵢ
    @test sol2[end, end] ≈ vₑ

    # Check various forms of interpolation with the forwards/backwards solution
    t = LinRange(sol1.t[1], sol1.t[2], 11)
    @test sol1(t[3], idxs=13) == sol2(t[3], idxs=13)
    @test sol1(t, idxs=13) == sol2(t, idxs=13)
    @test sol1(t[3], idxs=7:13) == sol2(t[3], idxs=7:13)
    @test sol1(t, idxs=7:13) == sol2(t, idxs=7:13)

    # Check that we can integrate orbital phase just as well
    sol3 = @test_logs min_level=Logging.Info inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Ω₁=Ωᵢ/2, Rᵢ=Rᵢ, integrate_orbital_phase=true, quiet=true)
    t₁, tₑ = extrema(sol3.t)
    t = sol2.t[t₁ .< sol2.t .< tₑ]
    @test sol2(t) ≈ sol3(t, idxs=1:13)
    @test sol3(0.0, idxs=14) ≈ 0.0  # Initial phase should be ≈0
    @test minimum(diff(sol3[end, :])) > 0  # Ensure that the phase is strictly increasing

end
