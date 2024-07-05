@testitem "Up-down instability" begin
    using Logging
    using Quaternionic

    # These parameters are taken more-or-less exactly from 064_CCE,
    # which is subject to this instability
    M₁ = 0.561844712025
    M₂ = 0.43822158103
    χ⃗₁ = QuatVec(-6.25908582173e-09, -2.21655853316e-08,  0.723734029888)
    χ⃗₂ = QuatVec( 3.39881185355e-08,  5.92671870568e-08, -0.799746224993)
    R = Rotor(1.0)
    Ωᵢ = 2e-3
    v = PostNewtonian.v(Ω=Ωᵢ, M=M₁+M₂)
    state = [M₁, M₂, vec(χ⃗₁)..., vec(χ⃗₂)..., components(R)..., v]

    Ω₊, Ω₋ = up_down_instability(state)

    # Test this general property, which should be true of *every* value returned for *any*
    # inputs to this function:
    @test 0 ≤ Ω₊ ≤ Ω₋ ≤ PostNewtonian.Ω(v=1, M=M₁+M₂)

    # Now, for the specific set of parameters given above, we should have
    @test Ω₊ ≈ 5.46e-4 rtol=1e-2
    @test Ω₋ ≈ PostNewtonian.Ω(v=1, M=M₁+M₂)

    ## Check on the test in `orbital_evolution` to make sure it works as expected
    # Ensure that even with `quiet=true`, we still get a warning for an unstable configuration
    udi_warning = r"\nThis system is likely to encounter the up-down instability in the"
    @test_logs (:warn, udi_warning) orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ω(v=v), quiet=true)

    # Ensure that incorrect mass ordering gets handled properly
    @test_logs (:warn, udi_warning) orbital_evolution(M₂, M₁, χ⃗₂, χ⃗₁, Ω(v=v), quiet=true)

    # Test that it doesn't warn if χ⃗₁ is downward or χ⃗₂ is upward
    @test_logs min_level=Logging.Info orbital_evolution(M₁, M₂, QuatVec(χ⃗₁.x, χ⃗₁.y, -χ⃗₁.z/2), χ⃗₂, Ω(v=v), Ωₑ=2Ω(v=v), quiet=true)
    @test_logs min_level=Logging.Info orbital_evolution(M₁, M₂, χ⃗₁, QuatVec(χ⃗₂.x, χ⃗₂.y, -χ⃗₂.z), Ω(v=v), quiet=true)

end
