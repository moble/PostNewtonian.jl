"""
    S⃗₁(pnsystem)

Dimensionful spin vector of object 1.
"""
S⃗₁(s::VecOrPNSystem) = χ⃗₁(s) * M₁(s)^2

"""
    S⃗₂(pnsystem)

Dimensionful spin vector of object 2.
"""
S⃗₂(s::VecOrPNSystem) = χ⃗₂(s) * M₂(s)^2

"""
    S⃗(pnsystem)
    S⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Total (dimensionful) spin vector ``S⃗₁+S⃗₂``.
"""
S⃗(M₁, M₂, χ⃗₁, χ⃗₂) = χ⃗₁ * M₁^2 + χ⃗₂ * M₂^2
S⃗(s::VecOrPNSystem) = S⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    Σ⃗(pnsystem)
    Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂)

Differential spin vector ``M(a⃗₂-a⃗₁)``.
"""
Σ⃗(M₁, M₂, χ⃗₁, χ⃗₂) =  (M₁ + M₂) * (χ⃗₂ * M₂ - χ⃗₁ * M₁)
Σ⃗(s::VecOrPNSystem) = Σ⃗(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗(pnsystem)
    χ⃗(S⃗, M)

Normalized spin vector ``S⃗/M²``.
"""
χ⃗(S⃗, M) = S⃗ / M^2
χ⃗(s::VecOrPNSystem) = χ⃗(S⃗(s), M(s))

"""
    χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂)

Symmetric spin vector ``(χ⃗₁+χ⃗₂)/2``.
"""
χ⃗ₛ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2
χ⃗ₛ(s::VecOrPNSystem) = χ⃗ₛ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

"""
    χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂)

Antisymmetric spin vector ``(χ⃗₁-χ⃗₂)/2``.
"""
χ⃗ₐ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
χ⃗ₐ(s::VecOrPNSystem) = χ⃗ₐ(M₁(s), M₂(s), χ⃗₁(s), χ⃗₂(s))

χₚₑᵣₚ(s::VecOrPNSystem) = √(χ₁²(s) - (χ₁ₗ(s))^2 + χ₂²(s) - (χ₂ₗ(s))^2)

χ₁²(s::VecOrPNSystem) = abs2vec(χ⃗₁(s))
χ₂²(s::VecOrPNSystem) = abs2vec(χ⃗₂(s))
χ₁(s::VecOrPNSystem) = absvec(χ⃗₁(s))
χ₂(s::VecOrPNSystem) = absvec(χ⃗₂(s))
χ₁₂(s::VecOrPNSystem) = χ⃗₁(s) ⋅ χ⃗₂(s)
χ₁ₗ(s::VecOrPNSystem) = χ⃗₁(s) ⋅ ℓ̂(s)
χ₂ₗ(s::VecOrPNSystem) = χ⃗₂(s) ⋅ ℓ̂(s)
χₛₗ(s::VecOrPNSystem) = χ⃗ₛ(s) ⋅ ℓ̂(s)
χₐₗ(s::VecOrPNSystem) = χ⃗ₐ(s) ⋅ ℓ̂(s)

Sₙ(s::VecOrPNSystem) = S⃗(s) ⋅ n̂(s)
Σₙ(s::VecOrPNSystem) = Σ⃗(s) ⋅ n̂(s)
Sλ(s::VecOrPNSystem) = S⃗(s) ⋅ λ̂(s)
Σλ(s::VecOrPNSystem) = Σ⃗(s) ⋅ λ̂(s)
Sₗ(s::VecOrPNSystem) = S⃗(s) ⋅ ℓ̂(s)
Σₗ(s::VecOrPNSystem) = Σ⃗(s) ⋅ ℓ̂(s)

S₁ₙ(s::VecOrPNSystem) = S⃗₁(s) ⋅ n̂(s)
S₁λ(s::VecOrPNSystem) = S⃗₁(s) ⋅ λ̂(s)
S₁ₗ(s::VecOrPNSystem) = S⃗₁(s) ⋅ ℓ̂(s)
S₂ₙ(s::VecOrPNSystem) = S⃗₂(s) ⋅ n̂(s)
S₂λ(s::VecOrPNSystem) = S⃗₂(s) ⋅ λ̂(s)
S₂ₗ(s::VecOrPNSystem) = S⃗₂(s) ⋅ ℓ̂(s)
