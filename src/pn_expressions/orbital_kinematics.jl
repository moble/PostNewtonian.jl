@doc raw"""
    γₚₙ(pnsystem)

Compute the post-Newtonian parameter
```math
\gamma_{\mathrm{PN}} \equiv \frac{G\, M}{r\, c^2},
```
where ``r`` is the magnitude of the orbital separation.  This quantity is counted as having
PN order 1, and is given by Eq. (4.3) of [Bohé et al.
(2013)](https://arxiv.org/abs/1212.5520) and Eq.  (3.32) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).  The latter only contains the spin-squared terms
for the non-precessing case.

Note that there is a 3PN gauge term of ``-22ν\\ln(r/r₀′)/3`` that is simply ignored here, as
it should cancel out of any physical quantity.
"""
@pn_expression function γₚₙ(pnsystem)
    v^2 * @pn_expansion(
        # Non-spinning terms; Eq. (4.3) of Bohé et al. (2013)
        1
        + v^2 * (1 - ν / 3)
        + v^4 * (1 - 65ν / 12)
        + v^6 * (1 + (-2203//2520 - 41π^2 / 192)ν + 229ν^2 / 36 + ν^3 / 81)

        # Spin-orbit terms; Eq. (4.3) of Bohé et al. (2013)
        + v^3 * (5//3 * sₗ + δ * σₗ)
        + v^5 * ((10//3 + 8ν/9) * sₗ + 2δ * σₗ)
        + v^7 * ((5 - 127ν/12 - 6ν^2) * sₗ + δ * (3 - 61ν/6 - 8ν^2/3) * σₗ)

        # Spin-squared terms; Eq. (3.32) of Bohé et al. (2015)
        + v^4 * (
            sₗ^2 * (-κ₊/2 - 1)
            + sₗ * σₗ * (-δ*κ₊/2 - δ + κ₋/2)
            + σₗ^2 * (δ*κ₋/4 - κ₊/4 + (κ₊/2 + 1)ν)
        )
        + v^6 * (
            sₗ^2 * (-11δ*κ₋/12 - 11κ₊/12 + 14//9 + (-κ₊/6 - 1//3)ν)
            + sₗ * σₗ * (5δ/3 + (-δ*κ₊/6 - δ/3 + 23κ₋/6)ν)
            + σₗ^2 * (1 + (δ*κ₋ - κ₊ - 2)ν + (κ₊/6 + 1//3)ν^2)
        )
    )
end


"""
    r(pnsystem)

Orbital separation given by the *coordinate* distance between the centers of objects in a
binary.

This is related to the post-Newtonian parameter [`γₚₙ`](@ref) by
```math
r = \frac{G\, M}{\gamma_{\mathrm{PN}}\, c^2}.
```
This value is the orbit-averaged quantity.
"""
r(pnsystem) = M(pnsystem) / γₚₙ(pnsystem)


"""
    ṙ(pnsystem)

Time derivative of the orbital separation.

This is an approximate value based only on the interactions between the spins of the binary
components, which would lead to ``ṙ = 0`` for non-precessing systems.  This value ignores
inspiral due to radiation reaction, and assumes that the precession timescale is much longer
than the orbital period.  It is given by Eq. (4.12a) of [Buonanno, Faye, Hinderer
(2013)](https://arxiv.org/abs/1209.6349)

See [`r`](@ref) for the description of the (orbit-averaged) orbital separation.
"""
@pn_expression function ṙ(pnsystem)
    - Ω / (2M^2 * r(pnsystem)) * ((n̂ ⋅ S⃗₀⁺) * (λ̂ ⋅ S⃗₀⁻) + (λ̂ ⋅ S⃗₀⁺) * (n̂ ⋅ S⃗₀⁻))
end
