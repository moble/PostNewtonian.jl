@doc raw"""
    TaylorT1!(u̇, p)

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
@pn_expression 2 function TaylorT1!(u̇, p)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(p)
    v̇ = - (𝓕(p) + Ṁ₁ + Ṁ₂) / 𝓔′(p)  # This expression is what makes this TaylorT1
    Ω⃗ = Ω⃗ₚ(p) + Ω * ℓ̂
    χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
    χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
    u̇[1] = Ṁ₁
    u̇[2] = Ṁ₂
    u̇[3:5] = vec((Ṡ₁ / M₁^2 - 2χ₁ * Ṁ₁/M₁) * χ̂₁ + Ω⃗ᵪ₁(p) × χ⃗₁)
    u̇[6:8] = vec((Ṡ₂ / M₂^2 - 2χ₂ * Ṁ₂/M₂) * χ̂₂ + Ω⃗ᵪ₂(p) × χ⃗₂)
    u̇[9:12] = components(Ω⃗ * R / 2)
    u̇[13] = v̇
    if length(u̇) > 13
        u̇[14] = Ω
    end
    nothing
end


@doc raw"""
    TaylorT4!(u̇, u, p, t)

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
"""
@pn_expression 3 function TaylorT4!(u̇, u, p, t)
    error("Not yet implemented")
end


@doc raw"""
    TaylorT5!(u̇, u, p, t)

Compute the right-hand side for the orbital evolution of a non-eccentric binary in the
"TaylorT5" approximant.

In this approximant, we compute ``\dot{v}`` by expanding the right-hand side of *the inverse
of* the usual expression
```math
\frac{1}{\dot{v}}=\frac{dt}{dv} = -\frac{\mathcal{E}'} {\mathcal{F} + \dot{M}_1 + \dot{M}_2}
```
as a series in ``v``, truncating again at the specified PN order, evaluating the result, and
then taking the inverse.  This approximant was introduced by [Ajith
(2011)](https://arxiv.org/abs/1107.1267).  Compare [`TaylorT1!`](@ref) and
[`TaylorT5!`](@ref).

Here, `u` is the ODE state vector, which should just refer to the `state` vector stored in
the [`PNSystem`](@ref) object `p`.  The parameter `t` represents the time, and will surely
always be unused in this package, but is part of the `DifferentialEquations` API.
"""
@pn_expression 3 function TaylorT5!(u̇, u, p, t)
    error("Not yet implemented")
end
