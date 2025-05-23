@doc raw"""
    γₚₙ₀(pnsystem, [r₀′])

This is a helper function for [`γₚₙ`](@ref); this computes the portion of `γₚₙ` that does
not depend on the logarithm of the separation.  Since `γₚₙ` is defined in terms of the
separation, that term makes the standard definition into an implicit function.

Specifically, `γₚₙ = γₚₙ₀ - (v/c)^8 * (22ln(r*c^2/(G*M)) / 3)ν`, where the argument of the
logarithm is precisely `1/γₚₙ` — by definition.  Note that there is a term just like the
second term in the expression that involves `r₀′`, which is a gauge parameter.  By default,
that term is simply ignored in this function, but if the optional argument is provided, it
is included.  This combination of the logarithms involving `r` and `r₀′` should drop out of
the result for any physical quantity.

"""
@pn_expression function γₚₙ₀(pnsystem, r₀′=0)
    lnr₀′c²╱GM = if iszero(r₀′)
        0
    else
        ln(r₀′ * c^2 / (G * M))
    end
    γ₀ = (v / c)^2 * @pn_expansion(
        # Non-spinning terms; Eq. (4.3) of Bohé et al. (2013)
        1 +
            (v / c)^2 * (1 - ν / 3) +
            (v / c)^4 * (1 - 65ν / 12) +
            (v / c)^6 * (
                1
                + (-2203//2520 - 41π^2 / 192 + 22lnr₀′c²╱GM / 3)ν
                + 229ν^2 / 36
                + ν^3 / 81
            )

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

            # Spin-cubed terms; Eq. (6.15) of Marsat (2015)
            +
            (v / c)^7 * (
                sₗ^3 * (-3κ₊ / 2 + λ₊ - 9) +
                sₗ^2 * σₗ * (-5δ * κ₊ / 2 + 3δ * λ₊ / 2 - 14δ + 3κ₋ - 3λ₋ / 2) +
                sₗ * σₗ^2 * (
                    13δ * κ₋ / 4 - 3δ * λ₋ / 2 - 13κ₊ / 4 + 3λ₊ / 2 - 5
                    + (11κ₊ / 2 - 3λ₊ + 29)ν
                ) +
                σₗ^3 * (
                    -5δ * κ₊ / 4 + δ * λ₊ / 2 + 5κ₋ / 4 - λ₋ / 2
                    + (δ * κ₊ - δ * λ₊ / 2 + 5δ - 7κ₋ / 2 + 3λ₋ / 2)ν
                )
            )
    )
end


@doc raw"""
    γₚₙ(pnsystem, [r₀′])
    inverse_separation(pnsystem, [r₀′])

Compute the post-Newtonian parameter
```math
\gamma_{\mathrm{PN}} \equiv \frac{G\, M}{r\, c^2},
```
where ``r`` is the magnitude of the orbital separation.  This quantity has PN order 1, and
is given by Eq. (4.3) of [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520), with
spin-squared terms from Eq.  (3.32) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529) and spin-cubed terms from [Marsat
(2014)](https://arxiv.org/abs/1411.4118).

Note that there is a 3PN gauge term of ``-22ν\ln(r/r₀')/3``.  While this value should cancel
out of any physical quantity, it is included here for completeness.  Computing it requires a
few Newton steps to get the value of ``γ`` because the ``\ln(r)`` term depends on
``\gamma``.

Specifically, we use the helper function [`γₚₙ₀`](@ref) to write `γₚₙ = γₚₙ₀ + (v/c)^8 *
(22ln(γₚₙ)/3)ν`; given the value of `γₚₙ₀`, the purpose of this function is to determine
`γₚₙ`.

The default value of `r₀'` is `G*M/c^2`, which makes its contribution vanish.
"""
@pn_expression function γₚₙ(pnsystem, r₀′=0)
    γ₀ = γₚₙ₀(pnsystem, r₀′; pn_expansion_reducer)

    if pn_order(pnsystem) ≥ 3
        # Account for the 3PN gauge term.  Note that the coefficient of the logarithm is
        # too small for the Lambert W function to give us a useful result, so we just
        # do a few Newton steps to get the value of γ = γ₀ + (v/c)^8 * (22ln(γ) / 3)ν
        a = (v / c)^8 * (22ν / 3)
        Δγ = let γ₀ = sum(γ₀)
            Δγᵢ = zero(a*ln(γ₀))
            for i ∈ 1:20  # Limit the possible number of steps, just in case things break
                γᵢ = γ₀ + Δγᵢ
                δγ = - (Δγᵢ + a*ln(γᵢ)) / (1 + a / (γᵢ))
                Δγᵢ += δγ
                if abs(δγ) < 10eps(γᵢ)
                    break
                end
            end
            Δγᵢ
        end
        if isa(pn_expansion_reducer, Val{sum})
            return γ₀ + Δγ
        else
            # γ₀ will be a `PNExpansion`, and we need to add the `Δγ` to term 7
            coeffs = Tuple(
                i==7 ? c + Δγ : c
                for (i,c) ∈ enumerate(γ₀.coeffs)
            )
            return typeof(γ₀)(coeffs)
        end
    else
        return γ₀
    end
end
const inverse_separation = γₚₙ


@doc raw"""
    γₚₙ′(pnsystem)
    inverse_separation_deriv(pnsystem)

Compute the derivative of [`γₚₙ`](@ref) with respect to `v`.

Note that we ignore the `r₀′` term in this function; that constant is obviously independent
of `v`, though it is multiplied by `M`, which is not independent of `v`.  This dependence,
however, should be at a much higher PN order than is currently available in any case, so we
ignore it for simplicity.

This computation uses [`γₚₙ₀`](@ref) along with the following derivation:
```math
\begin{align*}
γₚₙ &= γₚₙ₀ + (v/c)^8 (22 \ln(γₚₙ)/3)ν \\
γₚₙ' &= γₚₙ₀' + 8(v/c)^7 (22 \ln(γₚₙ)/3)ν + (v/c)^8 (22 γₚₙ'/3γₚₙ)ν \\
γₚₙ' &= \frac{γₚₙ₀' + 8(v/c)^7 (22 \ln(γₚₙ)/3)ν} {1 - (v/c)^8 (22/3γₚₙ)ν}
\end{align*}
```

"""
@pn_expression function γₚₙ′(pnsystem)
    if !isa(pn_expansion_reducer, Val{sum})
        throw(ArgumentError(
            "`PostNewtonian.γₚₙ′` not implemented for `pn_expansion_reducer` types other"
            *" than `Val{sum}`.  (That is, you can't get individual terms out of this.)"
        ))
    end

    γ₀′ = γₚₙ₀′(pnsystem)

    if pn_order(pnsystem) ≥ 3
        let γₚₙ = γₚₙ(pnsystem)
            (γ₀′ + 8(v / c)^7 * (22ln(γₚₙ) / 3)ν) / (1 - (v / c)^8 * (22 / 3γₚₙ) * ν)
        end
    else
        γ₀′
    end
end
const inverse_separation_deriv = γₚₙ′


@doc raw"""
    γₚₙ₀′(pnsystem)

Helper function to compute [`γₚₙ′`](@ref).  This just computes the derivative of
[`γₚₙ₀`](@ref) with respect to `v`; `γₚₙ′` takes care of the extra complications arising
from the Newton iterations in [`γₚₙ′`](@ref).

"""
@generated function γₚₙ₀′(
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
    γₚₙ₀formula = γₚₙ₀(fdpnsystem; pn_expansion_reducer=Val(PNExpansionReducer))

    # Now we take the derivative of γₚₙ₀ with respect to v.
    γₚₙ₀′ = SVector(FastDifferentiation.derivative(γₚₙ₀formula, v(fdpnsystem)))

    # Turn that into an Expr (FD insists on making it a function)
    in_place = true
    init_with_zeros = false
    γₚₙ₀′expr = FastDifferentiation.make_Expr(γₚₙ₀′, vars, in_place, init_with_zeros)

    # Now, we use `MacroTools` to get the body of the function.
    γₚₙ₀′body = MacroTools.unblock(MacroTools.splitdef(γₚₙ₀′expr)[:body])

    # # At this point, the function is just a long series of statements inside an `@inbounds`
    # # block, which we will want later, but first we need to extract them.
    MacroTools.@capture(γₚₙ₀′body, @inbounds begin
        γₚₙ₀′statements__
    end) || throw(
        ArgumentError(
            "\nNo @inbounds block found in γₚₙ₀′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\nγₚₙ₀′($pnsystem)" *
            "\n```",
        ),
    )

    # The γₚₙ₀′statements are mostly what we want, except that the last line is a return
    # statement.  We want that result, but we don't to return it yet; we want to wrap that
    # result, so we just get that returned quantity here.
    MacroTools.@capture(γₚₙ₀′statements[end], return γₚₙ₀′return_) || throw(
        ArgumentError(
            "\nNo return statement found in γₚₙ′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\nγₚₙ₀′($pnsystem)" *
            "\n```",
        ),
    )
    γₚₙ₀′statements[end] = γₚₙ₀′return

    if PNExpansionReducer === identity
        # When `pn_expansion_reducer=Val(identity)` is passed, we return a PNExpansion
        NMax = Int(2PNOrder + 1)
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{$(length(γₚₙ₀′)),$(eltype(ST))}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(γₚₙ₀′statements...)
            end
            return PNExpansion{$(length(γₚₙ₀′)),$(eltype(ST)),$NMax}(Tuple(result))
        end
    else
        # Otherwise, FD produces a 1-tuple, so we just extract the value from that.
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{1,$(eltype(ST))}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(γₚₙ₀′statements...)
            end
            return result[1]
        end
    end
end


"""
    r(pnsystem, [r₀′])
    separation(pnsystem, [r₀′])

Compute the separation between the two black holes.  This is essentially the multiplicative
inverse of [`γₚₙ`](@ref), with some factors of `G` and `M` thrown in.
"""
@pn_expression function r(pnsystem, r₀′=0)
    let γₚₙ = γₚₙ(pnsystem, r₀′)
        return G * M / (γₚₙ * c^2)
    end
end
const separation = r


"""
    r′(pnsystem, [r₀′])
    separation_deriv(pnsystem, [r₀′])

Compute the derivative of the separation between the two black holes with respect to `v`.

Note that we ignore a derivative of `M` that appears in the `r₀′` term, as explained in
[`γₚₙ′`](@ref).
"""
@pn_expression function r′(pnsystem, r₀′=0)
    let γₚₙ = γₚₙ(pnsystem, r₀′), γₚₙ′ = γₚₙ′(pnsystem)
        -γₚₙ′ * G * M / (γₚₙ^2 * c^2)
    end
end
const separation_deriv = r′


"""
    ṙ(pnsystem, [r₀′])
    separation_dot(pnsystem, [r₀′])

Compute the derivative of the separation between the two black holes with respect to time.
"""
@pn_expression function ṙ(pnsystem, r₀′=0)
    let γₚₙ = γₚₙ(pnsystem, r₀′), γₚₙ′ = γₚₙ′(pnsystem), 𝓕 = 𝓕(pnsystem), 𝓔′ = 𝓔′(pnsystem)
        𝓕 * γₚₙ′ * G * M / (γₚₙ^2 * c^2 * 𝓔′)
    end
end
const separation_dot = ṙ


"""
    γₚₙ⁻¹(γ, pnsystem, [r₀′])
    inverse_separation_inverse(γ, pnsystem, [r₀′])

Return `v` such that `γₚₙ(pnsystem, r₀′) = γ` when `pnsystem` is evaluated at `v`.

Note that the value of `v` in the input `pnsystem` is ignored; you may use any value.  It
may also be convenient to know that you can set the value of `v` in `pnsystem` to the
returned value using `PostNewtonian.vindex` as in
```julia
pnsystem.state[PostNewtonian.vindex] = γₚₙ⁻¹(γ, pnsystem)
```
See also [`r⁻¹`](@ref).
"""
function γₚₙ⁻¹(γ, pnsystem, r₀′=0)
    if 2γ ≥ 1
        @error "Error in γₚₙ⁻¹" γ pnsystem
        throw(ArgumentError("γ=$γ ≥ 1/2 describes a binary that has already merged"))
    elseif γ ≤ 0
        @error "Error in γₚₙ⁻¹" γ pnsystem
        throw(ArgumentError("γ=$γ ≤ 0 is unphysical"))
    end

    pnsystemᵥ = deepcopy(pnsystem)

    function newton_step(v, γ, pnsystemᵥ)
        # We denote by pnystemᵥ the `pnsystem` with the value of `v` set to this function's
        # argument, the trial value `v`.  We take a Newton step to find the value of `v`
        # such that γₚₙ(pnsystemᵥ) = γ; that is, we're finding the root of
        #   f(v) = γₚₙ(pnsystemᵥ) - γ
        # with
        #   f′(v) = γₚₙ′(pnsystemᵥ)
        pnsystemᵥ.state[vindex] = v
        γᵥ = γₚₙ(pnsystemᵥ, r₀′)
        γᵥ′ = γₚₙ′(pnsystemᵥ)
        return -( (γᵥ - γ) / γᵥ′ )
    end

    # We can get an initial guess by solving the leading-order equation
    # γₚₙ = v² + v²^2 * (1 - ν / 3) for v.  The quadratic formula gives us
    #   v² = (-1 ± √(1 + 4(1 - ν / 3)γ)) / (2(1 - ν / 3))
    # We obviously want the result to be positive, so we take the positive root.  We don't
    # need to worry about signs in these square-roots because (1-ν/3) will always be
    # between 1 and 11/12, and γ will always be strictly positive.
    vᵢ = let ν = ν(pnsystem)
        √((-1 + √(1 + 4(1 - ν / 3)γ)) / (2(1 - ν / 3)))
    end

    # Now we just do a few Newton steps to get the value of v.
    maxNsteps = 50  # Limit the possible number of steps, just in case...
    for i ∈ 1:maxNsteps
        δvᵢ = newton_step(vᵢ, γ, pnsystemᵥ)
        vᵢ += δvᵢ
        if abs(δvᵢ) < 10eps(vᵢ)
            break
        end
        if i==maxNsteps
            @error "Failure in γₚₙ⁻¹: Failed to converge after $i iterations" γ pnsystem vᵢ
        end
    end

    return vᵢ
end
const inverse_separation_inverse = γₚₙ⁻¹


"""
    r⁻¹(r, pnsystem, [r₀′])
    separation_inverse(r, pnsystem, [r₀′])

Return `v` such that `r = r(v)` when `pnsystem` is evaluated at `v`.

Note that the value of `v` in the input `pnsystem` is ignored; you may use any value.  It
may also be convenient to know that you can set the value of `v` in `pnsystem` to the
returned value using `PostNewtonian.vindex` as in
```julia
pnsystem.state[PostNewtonian.vindex] = r⁻¹(r, pnsystem)
```
See also [`γₚₙ⁻¹`](@ref).
"""
function r⁻¹(r, pnsystem, r₀′=0)
    let c = 1, G = 1, M = M(pnsystem)
        γ = G * M / (r * c^2)
        v = γₚₙ⁻¹(γ, pnsystem, r₀′)
    end
end
const separation_inverse = r⁻¹


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

using PostNewtonian:
    @pn_expansion,
    @pn_expression,
    M,
    M₁,
    M₂,
    ν,
    δ,
    χ₁ₗ,
    χ₂ₗ,
    χ₁₂,
    Ω,
    type_converter,
    PNExpansionParameter

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
const separation = r

"""
    ṙ(pnsystem)

Eq. (4.12), computed as ṙ = (dE/dt) / (dE/dr), re-expanded and truncated.
"""
@pn_expression function ṙ(pnsystem)
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
    using PostNewtonian: @pn_expansion, @pn_expression, ṙ, Kidder1995

    rng = Random.Xoshiro(1234)
    for pnsystem ∈ (rand(rng, BBH) for _ ∈ 1:1_000)
        # We know that Larry's expression is outdated.  It may get more so as we include
        # newer PN terms, so this tolerance may need to be adjusted.  This is more of a
        # sanity check.
        @test Kidder1995.ṙ(pnsystem) ≈ ṙ(pnsystem) rtol = 0.03
    end
end

@testitem "separation_inverse" begin
    using Random
    using PostNewtonian: PostNewtonian, γₚₙ, γₚₙ⁻¹, M₁index, M₂index, v, r, r⁻¹

    rng = Random.Xoshiro(1234)
    for _ ∈ 1:100_000
        # First, create a random system.  Make it NSNS to ensure that as many code paths as
        # possible are tested.  Ensure that v≤1/2 to avoid cases where the system has
        # already merged.
        pnsystem = rand(rng, NSNS; v=rand(rng) / 2)

        # Test γ
        vᵧ = γₚₙ⁻¹(γₚₙ(pnsystem), pnsystem)
        @test abs(1 - vᵧ / v(pnsystem)) < 3eps(typeof(vᵧ))

        # Now perturb the masses just enough to ensure that the total mass is significantly
        # different from 1, but not so different as to mess with the tolerance.
        pnsystem.state[M₁index] *= 1.03
        pnsystem.state[M₂index] *= 1.09

        # And re-test with `r` instead of `γ`.
        vᵣ = r⁻¹(r(pnsystem), pnsystem)
        @test abs(1 - vᵣ / v(pnsystem)) < 3eps(typeof(vᵣ))
    end
end
