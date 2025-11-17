"""
    S⃗₁(pnsystem)

Dimensionful spin vector of object 1.
"""
@public S⃗₁(s::PNSystem) = χ⃗₁(s) * M₁(s)^2

"""
    S⃗₂(pnsystem)

Dimensionful spin vector of object 2.
"""
@public S⃗₂(s::PNSystem) = χ⃗₂(s) * M₂(s)^2

"""
    S⃗(pnsystem)
    S⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Total (dimensionful) spin vector ``S⃗₁+S⃗₂``.
"""
@public S⃗(M₁, M₂, χ⃗₁, χ⃗₂) = χ⃗₁ * M₁^2 + χ⃗₂ * M₂^2
S⃗(s::PNSystem) = S⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    Σ⃗(pnsystem)
    Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Differential spin vector ``M(a⃗₂-a⃗₁)``.
"""
@public Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂) = (M₁ + M₂) * (χ⃗₂ * M₂ - χ⃗₁ * M₁)
Σ⃗(s::PNSystem) = Σ⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗(pnsystem)
    χ⃗(S⃗, M)

Normalized spin vector ``S⃗/M²``.
"""
@public χ⃗(S⃗, M) = S⃗ / M^2
χ⃗(s::PNSystem) = χ⃗(S⃗(s), M(s))

"""
    χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂)

Symmetric spin vector ``(χ⃗₁+χ⃗₂)/2``.
"""
@public χ⃗ₛ(χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2
χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2
χ⃗ₛ(s::PNSystem) = χ⃗ₛ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂)

Antisymmetric spin vector ``(χ⃗₁-χ⃗₂)/2``.
"""
@public χ⃗ₐ(χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
χ⃗ₐ(s::PNSystem) = χ⃗ₐ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

@public χₚₑᵣₚ(s::PNSystem) = √(χ₁²(s) - (χ₁ₗ(s))^2 + χ₂²(s) - (χ₂ₗ(s))^2)
@public const chi_perp = χₚₑᵣₚ

@doc raw"""
    χₑ(s)
    chi_eff(s)

Effective spin parameter of the system.

Defined as
```math
\chi_{\mathrm{eff}}
\colonequals \frac{c}{G} \left(
    \frac{\mathbf{S}_1}{M_1} + \frac{\mathbf{S}_2}{M_2}
\right) \cdot \frac{\hat{\mathbf{L}}_{\mathrm{N}}} {M}.
```
"""
@public χₑ(s::PNSystem) = (S₁ₗ(s) / M₁(s) + S₂ₗ(s) / M₂(s)) / M(s)
@public const chi_eff = χₑ

@doc raw"""
    χₚ(s)
    chi_p(s)

Effective precession spin parameter of the system.

Note that there are *two different* definitions of this quantity found in the literature.
The [original definition](https://arxiv.org/abs/1408.1810) (converted to the convention
where ``M_1 \geq M_2``) is
```math
\begin{gathered}
A_1 = 2 + \frac{3M_2}{2M_1} \\
A_2 = 2 + \frac{3M_1}{2M_2} \\
\chi_{\mathrm{p}} \colonequals \frac{1}{A_1 M_1^2}
\mathrm{max}\left(A_1 S_{1\perp}, A_2 S_{2\perp} \right).
\end{gathered}
```
In a paper from early in the detection era, the [LIGO collaboration used this
definition](https://arxiv.org/abs/1606.01210).

However, a [more recent paper](https://arxiv.org/abs/2010.04131) redefines this essentially
as ``M_1^2`` times that quantity.  Using the convention that ``q = M_2/M_1 \leq 1``, the
definition may be more compactly written as
```math
\chi_{\mathrm{p}} \colonequals \mathrm{max} \left(
    \chi_{1\perp}, \chi_{2\perp} q \frac{4q+3}{4+3q}
\right).
```
Again, a more recent paper by [LIGO/Virgo/KAGRA](https://arxiv.org/abs/2111.03606) uses this
convention.

Because it seems to be the trend, this function uses the latter definition.
"""
@public function χₚ(s::PNSystem)
    χ₁ₚₑᵣₚ = √(χ₁ₙ(s)^2 + χ₁λ(s)^2)
    χ₂ₚₑᵣₚ = √(χ₂ₙ(s)^2 + χ₂λ(s)^2)
    let q = 1 / q(s)  # This is to convert to LVK's convention
        if q > 1
            q, χ₁ₚₑᵣₚ, χ₂ₚₑᵣₚ = 1 / q, χ₂ₚₑᵣₚ, χ₁ₚₑᵣₚ
        end
        max(χ₁ₚₑᵣₚ, χ₂ₚₑᵣₚ * q * (4q + 3) / (4 + 3q))
    end
    # let M₁=M₁(s), M₂=M₂(s)
    #     if M₁ < M₂
    #         M₁, M₂ = M₂, M₁
    #     end
    #     B₁ = 2 + 3M₂ / 2M₁
    #     B₂ = 2 + 3M₁ / 2M₂
    #     S₁ₚₑᵣₚ = √(S₁ₙ(s)^2 + S₁λ(s)^2)
    #     S₂ₚₑᵣₚ = √(S₂ₙ(s)^2 + S₂λ(s)^2)
    #     max(B₁*S₁ₚₑᵣₚ, B₂*S₂ₚₑᵣₚ) / (B₁*M₁^2)
    # end
end
const chi_p = χₚ

@doc raw"""
    S⃗₀⁺(s)
    S⃗₀⁺(M₁, M₂, κ₁, κ₂, S⃗₁, S⃗₂)

Defined in Eq. (3.4) of [Buonanno et al. (2012)](https://arxiv.org/abs/1209.6349):
```math
\vec{S}_0^+
= \frac{M}{M_1} \left( \frac{\kappa_1} {\kappa_2} \right)^{1/4}
  \left( 1 + \sqrt{1 - \kappa_1 \kappa_2} \right)^{1/2} \vec{S}_1
  + \frac{M}{M_2} \left( \frac{\kappa_2} {\kappa_1} \right)^{1/4}
    \left( 1 - \sqrt{1 - \kappa_1 \kappa_2} \right)^{1/2} \vec{S}_2.
```
Note that, currently, ``\kappa_1`` and ``\kappa_2`` are both assumed to be equal to 1, as is
the case for black holes.  You can define `κ₁` and `κ₂` to have other values for your own
`PNSystem` types, and this function will work appropriately.

See also [`S⃗₀⁻`](@ref).
"""
@public S⃗₀⁺(s::PNSystem) = S⃗₀⁺(M₁(s), M₂(s), κ₁(s), κ₂(s), S⃗₁(s), S⃗₂(s))
function S⃗₀⁺(M₁, M₂, κ₁, κ₂, S⃗₁, S⃗₂)
    M = M₁ + M₂
    κᵣ = (κ₁ / κ₂)^(1//4)
    return (M / M₁) * κᵣ * √(1 + √(1 - κ₁ * κ₂)) * S⃗₁ +
           (M / M₂) / κᵣ * √(1 - √(1 - κ₁ * κ₂)) * S⃗₂
end
@public S₀⁺ₙ(s::PNSystem) = S⃗₀⁺(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁ₙ(s), S₂ₙ(s))
@public S₀⁺λ(s::PNSystem) = S⃗₀⁺(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁λ(s), S₂λ(s))
@public S₀⁺ₗ(s::PNSystem) = S⃗₀⁺(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁ₗ(s), S₂ₗ(s))

@doc raw"""
    S⃗₀⁻(s)
    S⃗₀⁻(M₁, M₂, κ₁, κ₂, S⃗₁, S⃗₂)

Defined below Eq. (3.4) of [Buonanno et al. (2012)](https://arxiv.org/abs/1209.6349):
```math
\vec{S}_0^-
= \frac{M}{M_1} \left( \frac{\kappa_1} {\kappa_2} \right)^{1/4}
  \left( 1 - \sqrt{1 - \kappa_1 \kappa_2} \right)^{1/2} \vec{S}_1
  + \frac{M}{M_2} \left( \frac{\kappa_2} {\kappa_1} \right)^{1/4}
    \left( 1 + \sqrt{1 - \kappa_1 \kappa_2} \right)^{1/2} \vec{S}_2.
```
Note that, currently, ``\kappa_1`` and ``\kappa_2`` are both assumed to be equal to 1, as is
the case for black holes.  You can define `κ₁` and `κ₂` to have other values for your own
`PNSystem` types, and this function will work appropriately.

See also [`S⃗₀⁺`](@ref).
"""
@public S⃗₀⁻(s::PNSystem) = S⃗₀⁻(M₁(s), M₂(s), κ₁(s), κ₂(s), S⃗₁(s), S⃗₂(s))
S⃗₀⁻(M₁, M₂, κ₁, κ₂, S⃗₁, S⃗₂) = S⃗₀⁺(M₂, M₁, κ₂, κ₁, S⃗₂, S⃗₁)
@public S₀⁻ₙ(s::PNSystem) = S⃗₀⁻(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁ₙ(s), S₂ₙ(s))
@public S₀⁻λ(s::PNSystem) = S⃗₀⁻(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁λ(s), S₂λ(s))
@public S₀⁻ₗ(s::PNSystem) = S⃗₀⁻(M₁(s), M₂(s), κ₁(s), κ₂(s), S₁ₗ(s), S₂ₗ(s))

@public χ₁²(s::PNSystem) = abs2vec(χ⃗₁(s))
@public χ₂²(s::PNSystem) = abs2vec(χ⃗₂(s))
@public χ₁(s::PNSystem) = absvec(χ⃗₁(s))
@public χ₂(s::PNSystem) = absvec(χ⃗₂(s))
@public χ₁₂(s::PNSystem) = χ⃗₁(s) ⋅ χ⃗₂(s)
@public χₛₗ(s::PNSystem) = χ⃗ₛ(s) ⋅ ℓ̂(s)
@public χₐₗ(s::PNSystem) = χ⃗ₐ(s) ⋅ ℓ̂(s)
@public χ₁ₙ(s::PNSystem) = χ⃗₁(s) ⋅ n̂(s)
@public χ₁λ(s::PNSystem) = χ⃗₁(s) ⋅ λ̂(s)
@public χ₁ₗ(s::PNSystem) = χ⃗₁(s) ⋅ ℓ̂(s)
@public χ₂ₙ(s::PNSystem) = χ⃗₂(s) ⋅ n̂(s)
@public χ₂λ(s::PNSystem) = χ⃗₂(s) ⋅ λ̂(s)
@public χ₂ₗ(s::PNSystem) = χ⃗₂(s) ⋅ ℓ̂(s)

@public Sₙ(s::PNSystem) = S⃗(s) ⋅ n̂(s)
@public Σₙ(s::PNSystem) = Σ⃗(s) ⋅ n̂(s)
@public Sλ(s::PNSystem) = S⃗(s) ⋅ λ̂(s)
@public Σλ(s::PNSystem) = Σ⃗(s) ⋅ λ̂(s)
@public Sₗ(s::PNSystem) = S⃗(s) ⋅ ℓ̂(s)
@public Σₗ(s::PNSystem) = Σ⃗(s) ⋅ ℓ̂(s)
@public sₗ(s::PNSystem) = S⃗(s) ⋅ ℓ̂(s) / M(s)^2
@public σₗ(s::PNSystem) = Σ⃗(s) ⋅ ℓ̂(s) / M(s)^2

@public S₁ₙ(s::PNSystem) = S⃗₁(s) ⋅ n̂(s)
@public S₁λ(s::PNSystem) = S⃗₁(s) ⋅ λ̂(s)
@public S₁ₗ(s::PNSystem) = S⃗₁(s) ⋅ ℓ̂(s)
@public S₂ₙ(s::PNSystem) = S⃗₂(s) ⋅ n̂(s)
@public S₂λ(s::PNSystem) = S⃗₂(s) ⋅ λ̂(s)
@public S₂ₗ(s::PNSystem) = S⃗₂(s) ⋅ ℓ̂(s)
