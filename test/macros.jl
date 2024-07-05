let # ignore: # @testitem "macros" begin
    using Quaternionic
    using Symbolics
    using PostNewtonian: var_collect
    @test var_collect(:(v^2 * (10ν/3)), :v) == (2, :((0, 0, 1 * ((10ν) / 3))))
    @test var_collect(:(v^3 * x + 1), :v) == (3, :((1, 0, 0, 1x)))
    @test var_collect(:(v^3 * x - 1), :v) == (3, :((-(1), 0, 0, 1x)))
    @test var_collect(:(v^3/M * ϕ̇̂₁ - Ωₕ₁), :v) == (3, :((-Ωₕ₁, 0, 0, (1 / M) * ϕ̇̂₁)))
    @test var_collect(:(v ^ 4), :v) == (4, :((0, 0, 0, 0, 1)))
    @test var_collect(:(v ^ 4 * (1)), :v) == (4, :((0, 0, 0, 0, 1 * 1)))
    @test var_collect(:(v ^ 4 * -(1)), :v) == (4, :((0, 0, 0, 0, 1 * -(1))))
    @test var_collect(:(v ^ 4 / x), :v) == (4, :((0, 0, 0, 0, 1 / x)))
    @test var_collect(:(v ^ 4 * y / x), :v) == (4, :((0, 0, 0, 0, (1y) / x)))
    @test var_collect(:(v ^ 4 * -(y) / x), :v) == (4, :((0, 0, 0, 0, (1 * -y) / x)))
    @test var_collect(:(v ^ 4 * -(y * z) / x), :v) == (4, :((0, 0, 0, 0, (1 * -(y*z)) / x)))
    @test var_collect(:(a + v ^ 4), :v) == (4, :((a, 0, 0, 0, 1)))
    @test var_collect(:(a + v ^ 4 * (1)), :v) == (4, :((a, 0, 0, 0, 1 * 1)))
    @test var_collect(:(a + v ^ 4 * -(1)), :v) == (4, :((a, 0, 0, 0, 1 * -(1))))
    @test var_collect(:(a + v ^ 4 / x), :v) == (4, :((a, 0, 0, 0, 1 / x)))
    @test var_collect(:(a + v ^ 4 * y / x), :v) == (4, :((a, 0, 0, 0, (1y) / x)))
    @test var_collect(:(a + v ^ 4 * -(y) / x), :v) == (4, :((a, 0, 0, 0, (1 * -y) / x)))

    @test_broken var_collect(:(v^3/M * ϕ̇̂₁ - Ωₕ₁ + 1), :v) == (3, :((-Ωₕ₁+1, 0, 0, (1 / M) * ϕ̇̂₁)))

    @eval PostNewtonian PostNewtonian.@pn_expression function macro_test(pnstate)
        v*π√10
    end

    for T ∈ [Float32, Float64, BigFloat]
        M₁,M₂,χ⃗₁,χ⃗₂,R,v = T(1),T(1),QuatVec(T(0)),QuatVec(T(0)),Rotor(T(1)),T(1//17)
        pnstate = BBH(;M₁, M₂, χ⃗₁, χ⃗₂, R, v)
        @test PostNewtonian.macro_test(pnstate) == v*(T(π)*√T(10))
    end

    @test isequal(
        PostNewtonian.macro_test(symbolic_pnsystem),
        PostNewtonian.v(symbolic_pnsystem)
        *Num(SymbolicUtils.Term(PostNewtonian.hold, [π]))
        *√Num(SymbolicUtils.Term(PostNewtonian.hold, [10]))
    )
end
