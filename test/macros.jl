@testset verbose=true "macros" begin
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
        PostNewtonian.v(symbolic_pnsystem)*Num(π)*√Num(SymbolicUtils.Term(identity, [10]))
    )

end
