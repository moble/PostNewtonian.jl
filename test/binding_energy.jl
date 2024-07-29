#= #

The most stringent test of the macros comes in `𝓔′`, because we have a typical PN
expression in `𝓔` (which is already a reasonable test of the macros), then we evaluate it
symbolically, differentiate that symbolically, convert it back to a piece of code, apply
`@pn_expansion`, and then wrap it up in a function to which we apply `@pn_expression` again.

So this test does all that a little more manually and compares the results at each order.

# =#

@testitem "binding_energy" begin
    using DoubleFloats: Double64

    include("binding_energy_reference.jl")

    for PNOrder ∈ 0//2 : 1//2 : 13//2

        for T ∈ [Float32, Float64, Double64]
            v = T(1//100)
            numpn = rand(NSNS; v, PNOrder)
            ϵ = 100eps(PostNewtonian.μ(numpn) * v^2)
            @test 𝓔(numpn) ≈ be(numpn, false) atol=ϵ rtol=100eps(T)
            @test 𝓔′(numpn) ≈ be(numpn, true) atol=ϵ rtol=100eps(T)
        end

    end

end

@testitem "binding_energy_symbolics" begin
    using Symbolics
    using DoubleFloats: Double64
    using Random: Xoshiro

    include("binding_energy_reference.jl")

    rng = Xoshiro(1234)

    for PNOrder ∈ [0//2 : 1//2 : 13//2; 1_000//2]

        for T ∈ [Float32, Float64, Double64]
            v = T(1//10)
            for _ ∈ 1:100
                numpn = rand(rng, NSNS; v, PNOrder)
                ϵ = 2eps(PostNewtonian.μ(numpn) * v^2)
                @test 𝓔′(numpn, Val(:Symbolics)) ≈ 𝓔′(numpn) atol=ϵ rtol=3eps(T)
            end
        end

    end

end
