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


χ₁²(s::PNState) = abs2vec(χ⃗₁(s))
χ₂²(s::PNState) = abs2vec(χ⃗₂(s))
χ₁(s::PNState) = absvec(χ⃗₁(s))
χ₂(s::PNState) = absvec(χ⃗₂(s))
χ₁₂(s::PNState) = χ⃗₁(s) ⋅ χ⃗₂(s)
χ₁ₗ(s::PNState) = χ⃗₁(s) ⋅ ℓ̂(s)
χ₂ₗ(s::PNState) = χ⃗₂(s) ⋅ ℓ̂(s)
χₛₗ(s::PNState) = χ⃗ₛ(s) ⋅ ℓ̂(s)
χₐₗ(s::PNState) = χ⃗ₐ(s) ⋅ ℓ̂(s)

Sₙ(s::PNState) = S⃗(s) ⋅ n̂(s)
Σₙ(s::PNState) = Σ⃗(s) ⋅ n̂(s)
Sλ(s::PNState) = S⃗(s) ⋅ λ̂(s)
Σλ(s::PNState) = Σ⃗(s) ⋅ λ̂(s)
Sₗ(s::PNState) = S⃗(s) ⋅ ℓ̂(s)
Σₗ(s::PNState) = Σ⃗(s) ⋅ ℓ̂(s)

S₁ₙ(s::PNState) = S⃗₁(s) ⋅ n̂(s)
S₁λ(s::PNState) = S⃗₁(s) ⋅ λ̂(s)
S₁ₗ(s::PNState) = S⃗₁(s) ⋅ ℓ̂(s)
S₂ₙ(s::PNState) = S⃗₂(s) ⋅ n̂(s)
S₂λ(s::PNState) = S⃗₂(s) ⋅ λ̂(s)
S₂ₗ(s::PNState) = S⃗₂(s) ⋅ ℓ̂(s)
