S⃗₁(s::PNState) = χ⃗₁(s) * M₁(s)^2
S⃗₂(s::PNState) = χ⃗₂(s) * M₂(s)^2

"""
    S⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Total spin vector ``S⃗₁+S⃗₂``.
"""
S⃗(M₁, M₂, χ⃗₁, χ⃗₂) = χ⃗₁ * M₁^2 + χ⃗₂ * M₂^2
S⃗(s::PNState) = S⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Differential spin vector ``M(a⃗₂-a⃗₁)``.
"""
Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂) =  (M₁ + M₂) * (χ⃗₂ * M₂ - χ⃗₁ * M₁)
Σ⃗(s::PNState) = Σ⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗(S⃗, M)

Normalized spin vector ``S⃗/M²``.
"""
χ⃗(S⃗, M) = S⃗ / M^2
χ⃗(s::PNState) = χ⃗(S⃗(s), M(s))

"""
    χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂)

Symmetric spin vector ``(χ⃗₁+χ⃗₂)/2``.
"""
χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2
χ⃗ₛ(s::PNState) = χ⃗ₛ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂)

Antisymmetric spin vector ``(χ⃗₁-χ⃗₂)/2``.
"""
χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
χ⃗ₐ(s::PNState) = χ⃗ₐ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))


Sₙ(s) = S⃗(s) ⋅ n̂(s)
Σₙ(s) = Σ⃗(s) ⋅ n̂(s)
Sλ(s) = S⃗(s) ⋅ λ̂(s)
Σλ(s) = Σ⃗(s) ⋅ λ̂(s)
Sₗ(s) = S⃗(s) ⋅ ℓ̂(s)
Σₗ(s) = Σ⃗(s) ⋅ ℓ̂(s)
