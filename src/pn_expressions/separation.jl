@doc raw"""
    γₚₙ(pnsystem)
    inverse_separation(pnsystem)

Compute the post-Newtonian parameter
```math
\gamma_{\mathrm{PN}} \equiv \frac{G\, M}{r\, c^2},
```
where ``r`` is the magnitude of the orbital separation.  This quantity has PN order 1, and
is given by Eq. (4.3) of [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520) and Eq.
(3.32) of [Bohé et al.  (2015)](https://arxiv.org/abs/1501.01529).

Note that there is a 3PN gauge term of ``-22ν\\ln(r/r₀′)/3`` that is simply ignored here, as
it should cancel out of any physical quantity.
"""
@pn_expression function γₚₙ(pnsystem)
    return (v / c)^2 * @pn_expansion(
        # Non-spinning terms; Eq. (4.3) of Bohé et al. (2013)
        1 +
            (v / c)^2 * (1 - ν / 3) +
            (v / c)^4 * (1 - 65ν / 12) +
            (v / c)^6 * (1 + (-2203//2520 - 41π^2 / 192)ν + 229ν^2 / 36 + ν^3 / 81)

            # Spin-orbit terms; Eq. (4.3) of Bohé et al. (2013)
            +
            (v / c)^3 * (5//3 * sₗ + δ * σₗ) +
            (v / c)^5 * ((10//3 + 8ν / 9) * sₗ + 2δ * σₗ) +
            (v / c)^7 * ((5 - 127ν / 12 - 6ν^2) * sₗ + δ * (3 - 61ν / 6 - 8ν^2 / 3) * σₗ)

            # Spin-squared terms; Eq. (3.32) of Bohé et al. (2015)
            +
            (v / c)^4 * (
                sₗ^2 * (-κ₊ / 2 - 1) +
                sₗ * σₗ * (-δ * κ₊ / 2 - δ + κ₋ / 2) +
                σₗ^2 * (δ * κ₋ / 4 - κ₊ / 4 + (κ₊ / 2 + 1)ν)
            ) +
            (v / c)^6 * (
                sₗ^2 * (-11δ * κ₋ / 12 - 11κ₊ / 12 + 14//9 + (-κ₊ / 6 - 1//3)ν) +
                sₗ * σₗ * (5δ / 3 + (-δ * κ₊ / 6 - δ / 3 + 23κ₋ / 6)ν) +
                σₗ^2 * (1 + (δ * κ₋ - κ₊ - 2)ν + (κ₊ / 6 + 1//3)ν^2)
            )
    )
end
const inverse_separation = γₚₙ

@doc raw"""
    γₚₙ′(pnsystem)
    inverse_separation_deriv(pnsystem)

Compute the derivative of [`γₚₙ`](@ref) with respect to `v`.

"""
@generated function γₚₙ′(
    pnsystem::PNSystem{ST,PNOrder}; pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum)
) where {ST,PNOrder,PNExpansionReducer}
    # Create a `PNSystem` with `FastDifferentiation` (henceforth FD) variables, using the
    # same PNOrder as the input `pnsystem`.
    fdpnsystem = FDPNSystem(eltype(ST), PNOrder)

    # FD expects a single vector of variables, so we concatenate the state vector with the
    # two tidal-coupling parameters
    vars = FastDifferentiation.Node[fdpnsystem.state; Λ₁(fdpnsystem); Λ₂(fdpnsystem)]

    # Now we evaluate γₚₙ using the FD variables.  This will expand all derived variables in
    # terms of the fundamental variables, but FD will take care of evaluating those
    # efficiently via common subexpression elimination (CSE).
    γₚₙformula = γₚₙ(fdpnsystem; pn_expansion_reducer=Val(PNExpansionReducer))

    # Now we take the derivative of γₚₙ with respect to v.
    γₚₙ′ = SVector(FastDifferentiation.derivative(γₚₙformula, v(fdpnsystem)))

    # Turn that into an Expr (FD insists on making it a function)
    in_place = true
    init_with_zeros = false
    γₚₙ′expr = FastDifferentiation.make_Expr(γₚₙ′, vars, in_place, init_with_zeros)

    # Now, we use `MacroTools` to get the body of the function.
    γₚₙ′body = MacroTools.unblock(MacroTools.splitdef(γₚₙ′expr)[:body])

    # # At this point, the function is just a long series of statements inside an `@inbounds`
    # # block, which we will want later, but first we need to extract them.
    MacroTools.@capture(γₚₙ′body, @inbounds begin
        γₚₙ′statements__
    end) || throw(
        ArgumentError(
            "\nNo @inbounds block found in γₚₙ′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\nγₚₙ′($pnsystem)" *
            "\n```",
        ),
    )

    # The γₚₙ′statements are mostly what we want, except that the last line is a return
    # statement.  We want that result, but we don't to return it yet; we want to wrap that
    # result, so we just get that returned quantity here.
    MacroTools.@capture(γₚₙ′statements[end], return γₚₙ′return_) || throw(
        ArgumentError(
            "\nNo return statement found in γₚₙ′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\nγₚₙ′($pnsystem)" *
            "\n```",
        ),
    )
    γₚₙ′statements[end] = γₚₙ′return

    if PNExpansionReducer === identity
        # When `pn_expansion_reducer=Val(identity)` is passed, we return a PNExpansion
        NMax = Int(2PNOrder + 1)
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{$(length(γₚₙ′)),$(eltype(ST))}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(γₚₙ′statements...)
            end
            return PNExpansion{$(length(γₚₙ′)),$(eltype(ST)),$NMax}(Tuple(result))
        end
    else
        # Otherwise, FD produces a 1-tuple, so we just extract the value from that.
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{1,$(eltype(ST))}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(γₚₙ′statements...)
            end
            return result[1]
        end
    end
end
const inverse_separation_deriv = γₚₙ′

"""
    r(pnsystem)
    separation(pnsystem)

Compute the separation between the two black holes.  This is essentially the inverse of
[`γₚₙ`](@ref), with some factors of `G` and `M` thrown in.

Note that there should be a factor of `1/c^2` in this expression; we reserve it to use
explicitly in PN expansions.  That is, for every factor of `1/r`, we explicitly include a
factor of `1/c^2` in the expansion.
"""
@pn_expression function r(pnsystem)
    let γₚₙ = γₚₙ(pnsystem)
        return G * M / (γₚₙ)
    end
end
const separation = r

"""
    r′(pnsystem)
    separation_prime(pnsystem)

Compute the derivative of the separation between the two black holes with respect to `v`.
"""
@pn_expression function r′(pnsystem)
    let γₚₙ = γₚₙ(pnsystem)
        -γₚₙ′ * G * M / (γₚₙ^2 * c^2)
    end
end
const separation_prime = r′

@pn_expression function ṙ(pnsystem)
    let γₚₙ = γₚₙ(pnsystem), γₚₙ′ = γₚₙ′(pnsystem), 𝓕 = 𝓕(pnsystem), 𝓔′ = 𝓔′(pnsystem)
        𝓕 * γₚₙ′ * G * M / (γₚₙ^2 * c^2 * 𝓔′)
    end
end
const separation_dot = ṙ

"""
This module contains a few expressions from [Kidder
(1995)](https://arxiv.org/abs/gr-qc/9506022).

This is mostly here for testing, because these expressions are not directly used in this
package: they are somewhat outdated and describe quantities that are not actually used in
this formulation.  However, they were used in the SpEC code as an initial guess for
eccentricity reduction, so we want to make sure that results from this package are
consistent with those from SpEC.

"""
module Kidder1995

using PostNewtonian: @pn_expansion, @pn_expression, M, M₁, M₂, ν, δ, χ₁ₗ, χ₂ₗ, χ₁₂, Ω

"""
    r(pnsystem)

Eq. (4.13).
"""
@pn_expression function r(pnsystem)
    let m = M, m₁ = M₁, m₂ = M₂, η = ν, δm = δ * M, χ₁L̂ₙŝ₁ = χ₁ₗ, χ₂L̂ₙŝ₂ = χ₂ₗ, ω = Ω
        m *
        (m * ω)^(-2//3) *
        @pn_expansion(
            1 - 1//3 * (3 - η) * (m * ω)^(2//3) / c^2 -
            (
                1//3 *
                ((χ₁L̂ₙŝ₁ * (2 * m₁^2 / m^2 + 3η)) + (χ₂L̂ₙŝ₂ * (2 * m₂^2 / m^2 + 3η)))
            ) * (m * ω) / c^3 +
                (η * (19//4 + η / 9) - 1//2 * η * (χ₁₂ - 3χ₁L̂ₙŝ₁ * χ₂L̂ₙŝ₂)) *
            (m * ω)^(4//3) / c^4
        )
    end
end

"""
    drdt(pnsystem)

Eq. (4.12), computed as ṙ = (dE/dt) / (dE/dr), re-expanded and truncated.
"""
@pn_expression function drdt(pnsystem)
    let r = r(pnsystem)
        let m = M, m₁ = M₁, m₂ = M₂, η = ν, δm = δ * M, χ₁L̂ₙŝ₁ = χ₁ₗ, χ₂L̂ₙŝ₂ = χ₂ₗ
            -64//5 *
            η *
            (m / r)^3 *
            @pn_expansion(
                1 - 1//336 * (1751 + 588η) * (m / r) / c^2 -
                (
                    7//12 * (
                        (χ₁L̂ₙŝ₁ * (19 * m₁^2 / m^2 + 15η)) +
                        (χ₂L̂ₙŝ₂ * (19 * m₂^2 / m^2 + 15η))
                    ) - 4π
                ) * (m / r)^(3//2) / c^3 -
                    5//48 * η * (59χ₁₂ - 173χ₁L̂ₙŝ₁ * χ₂L̂ₙŝ₂) * (m / r)^2 / c^4
            )
        end
    end
end

end  # module Kidder1995

@testitem "separation" begin
    using Random
    using PostNewtonian: @pn_expansion, @pn_expression, separation_dot
    using PostNewtonian.Kidder1995: drdt

    rng = Random.Xoshiro(1234)
    for pnsystem ∈ (rand(rng, BBH) for _ ∈ 1:1_000)
        # We know that Larry's expression is old and inaccurate.  It may get more so as we
        # include newer PN terms, so this tolerance may need to be adjusted.  This is more
        # of a sanity check.
        @test drdt(pnsystem) ≈ separation_dot(pnsystem) rtol = 0.03
    end
end
