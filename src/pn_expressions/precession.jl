"""
    Ω⃗ₚ(pnsystem)

Compute the angular velocity of orbital precession.

This is the angular velocity *of* the orbital angular velocity direction unit vector ``ℓ̂``;
the time derivative of that *unit* vector is ``Ω⃗ₚ × ℓ̂``.

At the moment, this is computed solely by expressions from [Bohé et al.
(2013)](https://arxiv.org/abs/1212.5520).  See [`𝛡`](@ref) for details.

See also [`R`](@ref PostNewtonian.R).
"""
Ω⃗ₚ(pnsystem) = 𝛡(pnsystem)
const Omega_p = Ω⃗ₚ


"""
    𝛡(pnsystem)

Compute the angular velocity of orbital precession according to Bohé et al.

As [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520) explain above their Eq. (4.1), the
orbital precession is given by the time derivative of the orbital axis: 𝓵̇ = 𝛡 × 𝓵, where
the angular velocity is along the separation vector 𝓷, so that 𝛡 = ϖ 𝓷.  And in turn,
they define aₗ ≔ r ω ϖ, where r is the separation and ω is the orbital angular frequency.
Then, they define the PN parameter γₚₙ≔M/r and we have Mω = v³ so that ϖ = γₚₙ aₗ / v³.  The
parameters γₚₙ and aₗ are given by Eqs. (4.3) and (4.4), and given here by the functions
[`γₚₙ`](@ref) and [`aₗ`](@ref).

The spin-squared terms (by which we mean both spin-spin and spin-orbit squared terms) in the
energy are known to 3pN order, and given in [Eq. (3.32) of Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).

See also [`R`](@ref PostNewtonian.R).
"""
@pn_expression function 𝛡(pnsystem)
    (γₚₙ(pnsystem) * aₗ(pnsystem) / (v/c)^3) * n̂
end


@doc raw"""
    γₚₙ(pnsystem)

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
    (v/c)^2 * @pn_expansion(
        # Non-spinning terms; Eq. (4.3) of Bohé et al. (2013)
        1
        + (v/c)^2 * (1 - ν / 3)
        + (v/c)^4 * (1 - 65ν / 12)
        + (v/c)^6 * (1 + (-2203//2520 - 41π^2 / 192)ν + 229ν^2 / 36 + ν^3 / 81)

        # Spin-orbit terms; Eq. (4.3) of Bohé et al. (2013)
        + (v/c)^3 * (5//3 * sₗ + δ * σₗ)
        + (v/c)^5 * ((10//3 + 8ν/9) * sₗ + 2δ * σₗ)
        + (v/c)^7 * ((5 - 127ν/12 - 6ν^2) * sₗ + δ * (3 - 61ν/6 - 8ν^2/3) * σₗ)

        # Spin-squared terms; Eq. (3.32) of Bohé et al. (2015)
        + (v/c)^4 * (
            sₗ^2 * (-κ₊/2 - 1)
            + sₗ * σₗ * (-δ*κ₊/2 - δ + κ₋/2)
            + σₗ^2 * (δ*κ₋/4 - κ₊/4 + (κ₊/2 + 1)ν)
        )
        + (v/c)^6 * (
            sₗ^2 * (-11δ*κ₋/12 - 11κ₊/12 + 14//9 + (-κ₊/6 - 1//3)ν)
            + sₗ * σₗ * (5δ/3 + (-δ*κ₊/6 - δ/3 + 23κ₋/6)ν)
            + σₗ^2 * (1 + (δ*κ₋ - κ₊ - 2)ν + (κ₊/6 + 1//3)ν^2)
        )
    )
end


"""
    aₗ(pnsystem)

Eq. (4.4) of [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520).  This term contributes
to [`𝛡`](@ref).
"""
@pn_expression function aₗ(pnsystem)
    (v/c)^7/M^3 * @pn_expansion(
        (7Sₙ + 3δ*Σₙ)
        + (v/c)^2 * ((-10 - 29ν/3) * Sₙ + δ*(-6 - 9ν/2) * Σₙ)
        + (v/c)^4 * (
            (3//2 + 59ν/4 + 52ν^2/9) * Sₙ
            + δ*(3//2 + 73ν/8 + 17ν^2/6) * Σₙ
        )
    )
end


"""
    Ω⃗ᵪ₁(pnsystem)

Compute the angular velocity of precession of χ⃗₁

In the approximation that the spin *magnitude* is constant, the time derivative of ``χ⃗₁``
is due to its rotation alone, and is given by ``Ω⃗ᵪ₁ × χ⃗₁``.

Note that this function simply calls [`Ω⃗ᵪ`](@ref) with the appropriate parameters.
"""
@pn_expression function Ω⃗ᵪ₁(pnsystem)
    QuatVec(Ω⃗ᵪ(M₁, M₂, χ⃗₁, χ⃗₂, v, R, pnsystem))
end
const Omega_chi1 = Ω⃗ᵪ₁


"""
    Ω⃗ᵪ₂(pnsystem)

Compute the angular velocity of precession of χ⃗₂

In the approximation that the spin *magnitude* is constant, the time derivative of ``χ⃗₂``
is due to its rotation alone, and is given by ``Ω⃗ᵪ₂ × χ⃗₂``.

Note that this function simply calls [`Ω⃗ᵪ`](@ref) with the appropriate parameters.
"""
@pn_expression function Ω⃗ᵪ₂(pnsystem)
    QuatVec(Ω⃗ᵪ(M₂, M₁, χ⃗₂, χ⃗₁, v, R, pnsystem))
end
const Omega_chi2 = Ω⃗ᵪ₂


"""
    Ω⃗ᵪ(Mⱼ, Mₖ, χ⃗ⱼ, χ⃗ₖ, R)

Compute the angular velocity of precession of spin vector `χ⃗ⱼ`.

In the approximation that the spin *magnitude* is constant, the time derivative of ``χ⃗ⱼ``
is due to its rotation alone, and is given by ``Ω⃗ᵪ × χ⃗ⱼ``.

Note that this function is called by [`Ω⃗ᵪ₁`](@ref) and [`Ω⃗ᵪ₂`](@ref) with the appropriate
parameters; you probably want to use those instead of this one.

The spin-spin term is given by Eq. (2.4) of [Kidder
(1995)](http://link.aps.org/doi/10.1103/PhysRevD.52.821); the spin-orbit terms by Eq. (4.5)
of [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2); and the quadrupole-monopole term
by Eq. (2.7) [Racine (2008)](http://link.aps.org/doi/10.1103/PhysRevD.78.044021).
"""
function Ω⃗ᵪ(Mⱼ, Mₖ, χ⃗ⱼ, χ⃗ₖ, v, R, pnsystem)
    # Note that we don't use the `@pn_expression` macro here, because we're swapping
    # the order of certain arguments above, so we do it manually here, and don't just
    # call things like `δ(pnsystem)` because that would fail to swap.
    let M=M(Mⱼ, Mₖ), ν=ν(Mⱼ, Mₖ), δ=δ(Mⱼ, Mₖ), n̂=n̂(R), ℓ̂=ℓ̂(R), c = one(eltype(pnsystem))
        χⱼₙ = χ⃗ⱼ ⋅ n̂
        χₖₙ = χ⃗ₖ ⋅ n̂

        (v/c)^5/M * @pn_expansion pnsystem (
            # Spin-spin term from Eq. (2.4) of Kidder
            v/c * (Mₖ^2 / M^2) * (-χ⃗ₖ + 3χₖₙ * n̂)

            # Spin-orbit terms from Eq. (4.5) of Bohé et al.
            + (
                (3//4 + ν/2 - 3δ/4)
                + (v/c)^2 * (9//16 + 5ν/4 - ν^2/24 + δ*(-9//16 + 5ν/8))
                + (v/c)^4 * (27//32 + 3ν/16 - 105ν^2/32 - ν^3/48 + δ*(-27//32 + 39ν/8 - 5ν^2/32))
            ) * ℓ̂

            # Quadrupole-monopole term from Eq. (2.7) of Racine
            + v/c * (3ν * χⱼₙ * n̂)
        )
    end
end
