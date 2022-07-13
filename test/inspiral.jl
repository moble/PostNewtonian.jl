@testset verbose=true "Inspiral" begin
    Random.seed!(1234)
    T, stepper = Float64, AutoVern9(Rodas5())
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

    sol1 = @test_logs (:info,forwards_termination) inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Rᵢ=Rᵢ)
    sol2 = @test_logs (:info,forwards_termination) (:info,backwards_termination) inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Ω₁=Ωᵢ/2, Rᵢ=Rᵢ)

    @test sol1.retcode == :Terminated
    @test sol1[end, 1] == vᵢ
    @test sol1[1] ≈ uᵢ
    @test sol1[end, end] ≈ vₑ

    @test sol2.retcode == :Terminated
    @test sol2[end, 1] ≈ v₁
    iᵢ = argmin(abs.(sol2.t))  # Assuming uᵢ corresponds to t==0
    @test sol2[iᵢ] ≈ uᵢ
    @test sol2[end, end] ≈ vₑ

end
