function v̇_numerator(pnsystem; pn_expansion_reducer=Val(sum))
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(pnsystem; pn_expansion_reducer)
    return -(𝓕(pnsystem; pn_expansion_reducer) + Ṁ₁ + Ṁ₂)
end

function v̇_denominator(pnsystem; pn_expansion_reducer=Val(sum))
    return 𝓔′(pnsystem; pn_expansion_reducer)
end

function v̇_numerator_coeffs(pnsystem)
    return v̇_numerator(pnsystem; pn_expansion_reducer=Val(identity)).coeffs
end

function v̇_denominator_coeffs(pnsystem)
    return v̇_denominator(pnsystem; pn_expansion_reducer=Val(identity)).coeffs
end

TaylorT1_v̇(p) = v̇_numerator(p) / v̇_denominator(p)
TaylorT4_v̇(p) = truncated_series_ratio(v̇_numerator_coeffs(p), v̇_denominator_coeffs(p))
function TaylorT5_v̇(p)
    return inv(truncated_series_ratio(v̇_denominator_coeffs(p), v̇_numerator_coeffs(p)))
end

"""
    causes_domain_error!(u̇, p)

Ensure that these parameters correspond to a physically valid set of PN parameters.

If the parameters are not valid, this function should modify `u̇` to indicate that the
current step is invalid.  This is done by filling `u̇` with `NaN`s, which will be detected
by the ODE solver and cause it to try a different (smaller) step size.

Currently, the only check that is done is to test that these parameters result in a PN
parameter v>0.  In the future, this function may be expanded to include other checks, or it
may be specialized for specific `PNSystem` subtypes.
"""
function causes_domain_error!(u̇::ST, p::PNSystem{NT}) where {ST,NT}
    if !ismutabletype(ST)
        error("`causes_domain_error!` cannot modify input `u̇` because it is immutable")
    end
    if v(p) ≤ 0  # If this is expanded, document the change in the docstring.
        u̇ .= convert(eltype(NT), NaN)
        true
    else
        false
    end
end

@pn_expression function TaylorTn!(pnsystem, u̇, TaylorTn_v̇::V̇) where {V̇}
    # If these parameters result in v≤0, fill u̇ with NaNs so that `solve` will
    # know that this was a bad step and try again.
    causes_domain_error!(u̇, pnsystem) && return nothing

    v̇ = TaylorTn_v̇(pnsystem)

    Ω⃗ = Ω⃗ₚ(pnsystem) + Ω * ℓ̂
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(pnsystem)
    χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
    χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
    χ⃗̇₁ = (Ṡ₁ / M₁^2) * χ̂₁ - (2Ṁ₁ / M₁) * χ⃗₁ + Ω⃗ᵪ₁(pnsystem) × χ⃗₁
    χ⃗̇₂ = (Ṡ₂ / M₂^2) * χ̂₂ - (2Ṁ₂ / M₂) * χ⃗₂ + Ω⃗ᵪ₂(pnsystem) × χ⃗₂
    Ṙ = Ω⃗ * R / 2
    u̇[M₁index] = Ṁ₁
    u̇[M₂index] = Ṁ₂
    u̇[χ⃗₁ˣindex] = χ⃗̇₁.x
    u̇[χ⃗₁ʸindex] = χ⃗̇₁.y
    u̇[χ⃗₁ᶻindex] = χ⃗̇₁.z
    u̇[χ⃗₂ˣindex] = χ⃗̇₂.x
    u̇[χ⃗₂ʸindex] = χ⃗̇₂.y
    u̇[χ⃗₂ᶻindex] = χ⃗̇₂.z
    u̇[Rʷindex] = Ṙ.w
    u̇[Rˣindex] = Ṙ.x
    u̇[Rʸindex] = Ṙ.y
    u̇[Rᶻindex] = Ṙ.z
    u̇[vindex] = v̇
    u̇[Φindex] = Ω
    return nothing
end

sys = SymbolCache(collect(pnsystem_symbols), nothing, :t)

@doc raw"""
    TaylorT1!(u̇, pnsystem)

Compute the right-hand side for the orbital evolution of a non-eccentric binary in the
"TaylorT1" approximant.

This approximant is the simplest, in which the time derivative ``\dot{v}`` is given directly
by
```math
\dot{v} = -\frac{\mathcal{F} + \dot{M}_1 + \dot{M}_2} {\mathcal{E}'},
```
and the PN expression for each term on the right-hand side is evaluated numerically before
insertion directly in this expression.  Compare [`TaylorT4!`](@ref) and [`TaylorT5!`](@ref).

Here, `u̇` is the time-derivative of the state vector, which is stored in the
[`PNSystem`](@ref) object `p`.
"""
TaylorT1!(u̇, pnsystem) = TaylorTn!(pnsystem, u̇, TaylorT1_v̇)
TaylorT1!(u̇, u, p, t) = (p.state.=u; TaylorT1!(u̇, p))

"""
    TaylorT1RHS!

A `SciMLBase.ODEFunction` wrapper for [`TaylorT1!`](@ref), suitable for passing into
`OrdinaryDiffEq.solve`.
"""
const TaylorT1RHS! = ODEFunction{true,FullSpecialize}(TaylorT1!; sys)

@doc raw"""
    TaylorT4!(u̇, pnsystem)

Compute the right-hand side for the orbital evolution of a non-eccentric binary in the
"TaylorT4" approximant.

In this approximant, we compute ``\dot{v}`` by expanding the right-hand side of
```math
\dot{v} = -\frac{\mathcal{F} + \dot{M}_1 + \dot{M}_2} {\mathcal{E}'}
```
as a series in ``v``, truncating again at the specified PN order, and only then is the
result evaluated.  Compare [`TaylorT1!`](@ref) and [`TaylorT5!`](@ref).

Here, `u` is the ODE state vector, which should just refer to the `state` vector stored in
the [`PNSystem`](@ref) object `p`.  The parameter `t` represents the time, and will surely
always be unused in this package, but is part of the `DifferentialEquations` API.

!!! note "Truncation order vs. `PNOrder` vs. PN order"
    When expanding the fraction given above as a series in ``v``, the truncation order is
    not necessarily the value of `PNOrder` given in the input `p`.  Instead, it is the
    highest order of the series that is present in the numerator or denominator — which is
    what we would normally *call* the PN order of those expansions.  The `PNOrder` parameter
    is the highest order of the series that is *allowed* to be present in those expansions,
    so that if `PNOrder` is `typemax(Int)`, the series will be expanded to the highest order
    given in any of the PN expansions, but the expansion of the ratio will not go to
    infinite order.  This is the reason that `TaylorT4` and `TaylorT5` do not approach
    `TaylorT1` as `PNOrder` approaches `typemax(Int)`.
"""
TaylorT4!(u̇, pnsystem) = TaylorTn!(pnsystem, u̇, TaylorT4_v̇)
TaylorT4!(u̇, u, p, t) = (p.state.=u; TaylorT4!(u̇, p))

"""
    TaylorT4RHS!

A `SciMLBase.ODEFunction` wrapper for [`TaylorT4!`](@ref), suitable for passing into
`OrdinaryDiffEq.solve`.
"""
const TaylorT4RHS! = ODEFunction{true,FullSpecialize}(TaylorT4!; sys)

@doc raw"""
    TaylorT5!(u̇, pnsystem)

Compute the right-hand side for the orbital evolution of a non-eccentric binary in the
"TaylorT5" approximant.

In this approximant, we compute ``\dot{v}`` by expanding the right-hand side of *the
multiplicative inverse of* the usual expression
```math
\frac{1}{\dot{v}}=\frac{dt}{dv} = -\frac{\mathcal{E}'} {\mathcal{F} + \dot{M}_1 + \dot{M}_2}
```
as a series in ``v``, truncating again at the specified PN order, evaluating the result, and
then taking the multiplicative inverse again to find ``\dot{v}``.  This approximant was
introduced by [Ajith (2011)](https://arxiv.org/abs/1107.1267) [see Eq. (3.5)].  Compare
[`TaylorT1!`](@ref) and [`TaylorT4!`](@ref).

Here, `u` is the ODE state vector, which should just refer to the `state` vector stored in
the [`PNSystem`](@ref) object `p`.  The parameter `t` represents the time, and will surely
always be unused in this package, but is part of the `DifferentialEquations` API.
"""
TaylorT5!(u̇, pnsystem) = TaylorTn!(pnsystem, u̇, TaylorT5_v̇)
TaylorT5!(u̇, u, p, t) = (p.state.=u; TaylorT5!(u̇, p))

"""
    TaylorT5RHS!

A `SciMLBase.ODEFunction` wrapper for [`TaylorT5!`](@ref), suitable for passing into
`OrdinaryDiffEq.solve`.
"""
const TaylorT5RHS! = ODEFunction{true,FullSpecialize}(TaylorT5!; sys)
