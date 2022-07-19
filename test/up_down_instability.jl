@testset verbose=true "Up-down instability" begin
    M₁=0.561844712025
    M₂=0.43822158103
    χ⃗₁ = QuatVec(-6.25908582173e-09, -2.21655853316e-08,  0.723734029888)
    χ⃗₂ = QuatVec( 3.39881185355e-08,  5.92671870568e-08, -0.799746224993)
    R = Rotor(1.0)
    v = 0.25
    u = [M₁, M₂, χ⃗₁.vec..., χ⃗₂.vec..., R.components..., v]

    Ω₊, Ω₋ = up_down_instability(u)
    @test 0 ≤ Ω₊ ≤ Ω₋ ≤ 1  # General property
    @test Ω₊ < 6e-4
    @test Ω₋ == 1

    udi_warning = r"This system is likely to encounter the up-down instability in the\n"
    @test_logs (:warn, udi_warning) inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ω(v=v), quiet=true)

end
