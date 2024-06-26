function v̇_numerator(p)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(p)
    - (𝓕(p) + Ṁ₁ + Ṁ₂)
end

function v̇_denominator(p)
    𝓔′(p)
end

function v̇_numerator_coeffs(p)
    error("Not yet implemented")
end

function v̇_denominator_coeffs(p)
    error("Not yet implemented")
end


# This represents most of the lines of code needed to compute each approximant, but is
# identical in each.  The only things that this code don't do are (1) checking
# `causes_domain_error!` and (2) computing v̇.  The first has to be done before the second,
# but is just one line, and the second changes depending on the approximant.
RHS_body = quote
    Ω⃗ = Ω⃗ₚ(p) + Ω * ℓ̂
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(p)
    χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
    χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
    χ⃗̇₁ = (Ṡ₁ / M₁^2) * χ̂₁ - (2Ṁ₁ / M₁) * χ⃗₁ + Ω⃗ᵪ₁(p) × χ⃗₁
    χ⃗̇₂ = (Ṡ₂ / M₂^2) * χ̂₂ - (2Ṁ₂ / M₂) * χ⃗₂ + Ω⃗ᵪ₂(p) × χ⃗₂
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
    nothing
end

sys = SymbolCache(collect(pnsystem_symbols), nothing, :t)

@eval @doc raw"""
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
@pn_expression 2 function TaylorT1!(u̇, p)
    # If these parameters result in v≤0, fill u̇ with NaNs so that `solve` will
    # know that this was a bad step and try again.
    causes_domain_error!(u̇, p) && return

    # This expression is what makes this TaylorT1
    v̇ = v̇_numerator(p) / v̇_denominator(p)

    $RHS_body
end

const TaylorT1RHS! = ODEFunction{true, FullSpecialize}(
    (u̇,u,p,t) -> (p.state.=u; TaylorT1!(u̇,p)); sys
)


@eval @doc raw"""
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
"""
@pn_expression 2 function TaylorT4!(u̇, p)
    # If these parameters result in v≤0, fill u̇ with NaNs so that `solve` will
    # know that this was a bad step and try again.
    causes_domain_error!(u̇, p) && return

    # This expression is what makes this TaylorT4
    v̇ = truncated_series_ratio(v̇_numerator_coeffs(p), v̇_denominator_coeffs(p), v)

    $RHS_body
end

const TaylorT4RHS! = ODEFunction{true, FullSpecialize}(
    (u̇,u,p,t) -> (p.state.=u; TaylorT4!(u̇,p)); sys
)


@eval @doc raw"""
    TaylorT5!(u̇, pnsystem)

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
@pn_expression 2 function TaylorT5!(u̇, p)
    # If these parameters result in v≤0, fill u̇ with NaNs so that `solve` will
    # know that this was a bad step and try again.
    causes_domain_error!(u̇, p) && return

    # This expression is what makes this TaylorT5
    v̇ = inv(truncated_series_ratio(v̇_denominator_coeffs(p), v̇_numerator_coeffs(p), v))

    $RHS_body
end

const TaylorT5RHS! = ODEFunction{true, FullSpecialize}(
    (u̇,u,p,t) -> (p.state.=u; TaylorT5!(u̇,p)); sys
)
