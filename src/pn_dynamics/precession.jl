"""
    Ω⃗ₚ(pn)

Compute the angular velocity of orbital precession.

This is the angular velocity *of* the orbital angular velocity direction unit
vector ``ℓ̂``; the time derivative of that *unit* vector is ``Ω⃗ₚ × ℓ̂``.

At the moment, this is computed solely by expressions from [Bohé et
al. (2013)](https://arxiv.org/abs/1212.5520).  See [`𝛡`](@ref) for details.

"""
Ω⃗ₚ(pn) = 𝛡(pn)
const Omega_p = Ω⃗ₚ


"""
    Ω⃗ᵪ₁(pn)

Compute the angular velocity of precession of χ⃗₁

In the approximation that the spin *magnitude* is constant, the time derivative
of ``χ⃗₁`` is due to its rotation alone, and is given by ``Ω⃗ᵪ₁ × χ⃗₁``.

Note that this function simply calls [`Ω⃗ᵪ`](@ref) with the appropriate parameters.

"""
function Ω⃗ᵪ₁(pn)
    Ω⃗ᵪ(pn.M₁, pn.M₂, pn.χ⃗₁, pn.χ⃗₂, pn.v, pn.R)
end
const Omega_chi1 = Ω⃗ᵪ₁


"""
    Ω⃗ᵪ₂(pn)

Compute the angular velocity of precession of χ⃗₂

In the approximation that the spin *magnitude* is constant, the time derivative
of ``χ⃗₂`` is due to its rotation alone, and is given by ``Ω⃗ᵪ₂ × χ⃗₂``.

Note that this function simply calls [`Ω⃗ᵪ`](@ref) with the appropriate parameters.

"""
function Ω⃗ᵪ₂(pn)
    Ω⃗ᵪ(pn.M₂, pn.M₁, pn.χ⃗₂, pn.χ⃗₁, pn.v, pn.R)
end
const Omega_chi2 = Ω⃗ᵪ₂


"""
    Ω⃗ᵪ(Mⱼ, Mₖ, χ⃗ⱼ, χ⃗ₖ, R)

Compute the angular velocity of precession of spin vector `χ⃗ⱼ`.

In the approximation that the spin *magnitude* is constant, the time derivative
of ``χ⃗ⱼ`` is due to its rotation alone, and is given by ``Ω⃗ᵪ × χ⃗ⱼ``.

Note that this function is called by [`Ω⃗ᵪ₁`](@ref) and [`Ω⃗ᵪ₂`](@ref) with the
appropriate parameters; you may want to use those instead.

The spin-spin term is given by Eq. (2.4) of [Kidder
(1995)](http://link.aps.org/doi/10.1103/PhysRevD.52.821); the spin-orbit terms
by Eq. (4.5) of [Bohé et al. (2013)](http://arxiv.org/abs/1212.5520v2); and the
quadrupole-monopole term by Eq. (2.7) [Racine
(2008)](http://link.aps.org/doi/10.1103/PhysRevD.78.044021).

"""
function Ω⃗ᵪ(Mⱼ, Mₖ, χ⃗ⱼ, χ⃗ₖ, v, R)
    M = Mⱼ + Mₖ
    let ν=ν(Mⱼ, Mₖ), δ=δ(Mⱼ, Mₖ), n̂=n̂(R), ℓ̂=ℓ̂(R)
        χⱼₙ = χ⃗ⱼ ⋅ n̂
        χₖₙ = χ⃗ₖ ⋅ n̂

        v^5/M * (
            # Spin-spin term from Eq. (2.4) of Kidder
            v * (Mₖ^2 / M^2) * (-χ⃗ₖ + 3χₖₙ * n̂)

            # Spin-orbit terms from Eq. (4.5) of Bohé et al.
            + (
                (3//4 + ν/2 - 3δ/4)
                + v^2 * (9//16 + 5ν/4 - ν^2/24 + δ*(-9//16 + 5ν/8))
                + v^4 * (27//32 + 3ν/16 - 105ν^2/32 - ν^3/48 + δ*(-27//32 + 39ν/8 - 5ν^2/32))
            ) * ℓ̂

            # Quadrupole-monopole term from Eq. (2.7) of Racine
            + v * (3ν * χⱼₙ * n̂)
        )
    end
end


"""
    𝛡(pn)

Compute the angular velocity of orbital precession according to Bohé et al.

As [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520) explain above their
Eq. (4.1), the orbital precession is given by the time derivative of the
orbital axis: ``𝓵̇ = 𝛡 × 𝓵``, where the angular velocity is along the separation
vector ``𝓷``, so that ``𝛡 = ϖ 𝓷``.  And in turn, they define ``aₗ ≔ r ω ϖ``,
where ``r`` is the separation and ``ω`` is the orbital angular frequency.
Then, they define the PN parameter ``γ≔M/r`` and we have ``Mω = v³`` so that
``ϖ = γ aₗ / v³``.  The parameters ``γ`` and ``aₗ`` are given by Eqs. (4.3) and
(4.4), and given here by the functions [`γ`](@ref) and  [`aₗ`](@ref).

"""
function 𝛡(pn)
    let n̂ = n̂(pn.R), v=pn.v
        (γ(pn) * aₗ(pn) / v^3) * n̂
    end
end


"""
    γ(pn)

Eq. (4.3) of [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520).  This term
contributes to [`𝛡`](@ref).

Note that there is a 3PN term of ``-22ν\\ln(r/r₀′)/3`` that is simply ignored
here.

"""
function γ(pn)
    @unpack pn
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), ℓ̂=ℓ̂(R), π=oftype(v, π)
        Sₗ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
        Σₗ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
        v^2 * (
            1
            + v^2 * (1 - ν / 3)
            + v^4 * (1 - 65ν / 12)
            + v^6 * (1 + (-2203//2520 - 41π^2 / 192)ν + 229ν^2 / 36 + ν^3 / 81)
            + v^3 * ((5//3 * Sₗ + δ * Σₗ) / M^2)
            + v^5 * (((10//3 + 8ν/9) * Sₗ + 2δ * Σₗ) / M^2)
            + v^7 * (((5 - 127ν/12 - 6ν^2) * Sₗ + δ * (3 - 61ν/6 - 8ν^2/3) * Σₗ)/M^2)
        )
    end
end


"""
    aₗ(pn)

Eq. (4.4) of [Bohé et al. (2013)](https://arxiv.org/abs/1212.5520).  This term
contributes to [`𝛡`](@ref).

"""
function aₗ(pn)
    @unpack pn
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), n̂=n̂(R)
        Sₙ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ n̂
        Σₙ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ n̂
        v^7/M^3 * (
            (7Sₙ + 3δ*Σₙ)
            + v^2 * ((-10 - 29ν/3) * Sₙ + δ*(-6 - 9ν/2) * Σₙ)
            + v^4 * (
                (3//2 + 59ν/4 + 52ν^2/9) * Sₙ
                + δ*(3//2 + 73ν/8 + 17ν^2/6) * Σₙ
            )
        )
    end
end

