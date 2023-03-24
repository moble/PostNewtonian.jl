#= #

The most stringent test of the macros comes in `𝓔′`, because we have a typical PN
expression in `𝓔` (which is already a reasonable test of the macros), then we evaluate it
symbolically, differentiate that symbolically, convert it back to a piece of code, apply
`@pn_expansion`, and then wrap it up in a function to which we apply `@pn_expression` again.

So this test does all that a little more manually and compares the results at each order.

# =#

@testset verbose=true "binding_energy" begin

using PostNewtonian: a₆ᶜ¹, a₆₅ᶜ¹, a₇ˡⁿ¹, a₇ᶜ¹

function be(pnsystem, deriv)

    let M₁ = PostNewtonian.M₁(pnsystem), M₂ = PostNewtonian.M₂(pnsystem),
        v = PostNewtonian.v(pnsystem),
        λ₁ = PostNewtonian.λ₁(pnsystem), λ₂ = PostNewtonian.λ₂(pnsystem),
        M = PostNewtonian.M(pnsystem), Sₗ = PostNewtonian.Sₗ(pnsystem),
        δ = PostNewtonian.δ(pnsystem), μ = PostNewtonian.μ(pnsystem),
        ν = PostNewtonian.ν(pnsystem),
        Σₗ = PostNewtonian.Σₗ(pnsystem), χ₁² = PostNewtonian.χ₁²(pnsystem),
        χ₁₂ = PostNewtonian.χ₁₂(pnsystem), χ₂² = PostNewtonian.χ₂²(pnsystem),
        χₐₗ = PostNewtonian.χₐₗ(pnsystem), χₛₗ = PostNewtonian.χₛₗ(pnsystem),
        π = PostNewtonian.type_converter(pnsystem, π),
        ln2 = PostNewtonian.type_converter(pnsystem, PostNewtonian.ln2),
        ln3 = PostNewtonian.type_converter(pnsystem, PostNewtonian.ln3),
        γₑ = PostNewtonian.type_converter(pnsystem, PostNewtonian.γₑ),
        ln = (x->log(PostNewtonian.type_converter(pnsystem, x))),
        pn_order=PostNewtonian.pn_order(pnsystem)

        e = Dict()
        eˡⁿ = Dict()

        c = -1//2 * μ  # NOTE: Included v^2 factor inside sum for easier differentiation
        e[0] = 1
        e[2] = (-ν/12 - 3//4)
        e[4] = (-ν^2/24 + 19ν/8 - 27//8)
        e[6] =  (-35ν^3/5184 - 155ν^2/96 + (34445//576 - 205π^2/96)ν - 675//64)
        e[8] = (
            -3969//128 + 77ν^4/31104 + 301ν^3/1728 + (-498449//3456 + 3157π^2/576)ν^2
            + (-123671//5760 + 1792ln2/15 + 9037π^2/1536 + 896γₑ/15)ν
        )
        eˡⁿ[8] = (448ν/15)

        # Below are the incomplete terms from Eq. (74) of https://arxiv.org/abs/1312.2503v2
        e[10] = (
            -45927//512 + ν^5/512 + 55ν^4/512 + (-1353π^2/256 + 69423//512)ν^3
            + (-21337π^2/1024 + 3a₆ᶜ¹ - 896ln2/5 - 448γₑ/5 + 893429//2880)ν^2
            + (-228916843//115200 - 9976γₑ/35 + 729ln3/7 - 23672ln2/35 + 126779π^2/512)ν
        )
        eˡⁿ[10] = (-4988ν/35 - 656ν^2/5)
        e[11] = (10ν/3 * (13696π/525 + ν*a₆₅ᶜ¹))
        e[12] = (
            -264627//1024
            + 2717ν^6/6718464
            + 5159ν^5/248832
            + (272855π^2/124416 - 20543435//373248)ν^4
            + (
                1232γₑ/27 + 6634243π^2/110592
                - 11a₆ᶜ¹/2 - 71700787//51840  + 2464ln2/27
            )ν^3
            + (
                113176680983//14515200 + 18491π^4/2304
                + 246004ln2/105 + 112772γₑ/105 + a₆ᶜ¹*11//2 + a₇ˡⁿ¹*2//3
                + a₇ᶜ¹*11//3 - 86017789π^2/110592 - 2673ln3/14
            )ν^2
            + (
                -389727504721//43545600 + 74888ln2/243 - 7128ln3/7
                - 30809603π^4/786432 - 3934568γₑ/8505 + 9118627045π^2/5308416
            )ν
        )
        eˡⁿ[12] = (
            - 1967284ν/8505
            + 24464ν^3/135
            + (39754//105 + a₇ˡⁿ¹*11//3)ν^2
        )

        # Spin-orbit
        e[3] = ((14Sₗ/3 + 2δ * Σₗ) / M^2)
        e[5] = (((11-61ν/9) * Sₗ + δ*(3 - 10ν/3) * Σₗ) / M^2)
        e[7] = (((135//4 - 367ν/4 + 29ν^2/12) * Sₗ + δ*(27//4 - 39ν + 5ν^2/4) * Σₗ) / M^2)

        # Spin-squared
        e[4] += (
            (1 + δ - 2ν) * (χ₁² + χ₂²)/4 - 3*(χₐₗ^2 + χₛₗ^2)/2
            - δ*(χ₂²/2 + 3χₐₗ*χₛₗ) + (χ₁₂ + 6χₐₗ^2)ν
        )

        # NS tidal coupling
        e[10] += (-9*((M₁/M₂)λ₂ + (M₂/M₁)λ₁) / M^5)
        e[12] += (
            (
                -11//2*(M₁/M₂)*(3+2M₂/M+3*(M₂/M)^2)λ₂
                - 11//2*(M₂/M₁)*(3+2M₁/M+3*(M₁/M)^2)λ₁
            ) / M^5
        )

        if deriv
            c * (
                sum((k+2)*v^(k+1)*coeff for (k,coeff) ∈ e if k ≤ 2pn_order; init=0)
                + sum(v^(k+1)*coeff*2*((k+2)*log(v)+1) for (k,coeff) ∈ eˡⁿ if k ≤ 2pn_order; init=0)
            )
        else
            c * v^2 * (
                sum(v^(k)*coeff for (k,coeff) ∈ e if k ≤ 2pn_order; init=0)
                + sum(v^(k)*coeff*2log(v) for (k,coeff) ∈ eˡⁿ if k ≤ 2pn_order; init=0)
            )
        end
    end
end

for PNOrder ∈ 0//2:1//2:15//2
    sympn = SymbolicPNSystem(PNOrder)

    𝓔1 = 𝓔(sympn)
    𝓔2 = be(sympn, false)
    diff = simplify(𝓔1-𝓔2, expand=true)
    @test iszero(diff)

    𝓔′1 = 𝓔′(sympn)
    𝓔′2 = be(sympn, true)
    diff′ = simplify(𝓔′1-𝓔′2, expand=true)
    @show PNOrder 𝓔′1 𝓔′2 diff′
    println()
    @test iszero(diff′)

    for T ∈ [Float32, Float64, Double64, BigFloat]
        v = T(1//100)
        pn_system = randn(NSNS; v, PNOrder)
        ϵ = 4eps(PostNewtonian.μ(pn_system) * v^2)
        @test 𝓔(pn_system) ≈ be(pn_system, false) atol=ϵ
        𝓔′3 = 𝓔′(pn_system)
        𝓔′4 = be(pn_system, true)
        @test 𝓔′3 ≈ 𝓔′4 atol=ϵ
        #@test 𝓔′(pn_system) ≈ be(pn_system, true) atol=ϵ
    end

end

end
