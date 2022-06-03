"""
    χ⃗(S⃗, M)

Normalized spin vector ``S⃗/M²``.

"""
χ⃗(S⃗, M) = S⃗ / M^2


"""
    S(M₁, M₂, χ⃗₁, χ⃗₂)

Total spin vector ``S⃗₁+S⃗₂``.

"""
S(M₁, M₂, χ⃗₁, χ⃗₂) = χ⃗₁ * M₁^2 + χ⃗₂ * M₂^2


"""
    Σ(M₁, M₂, χ⃗₁, χ⃗₂)

Differential spin vector ``M(a⃗₂-a⃗₁)``.

"""
Σ(M₁, M₂, χ⃗₁, χ⃗₂) =  (M₁ + M₂) * (χ⃗₂ * M₂ - χ⃗₁ * M₁)


"""
    χₛ(M₁, M₂, χ⃗₁, χ⃗₂)

Symmetric spin vector (χ⃗₁+χ⃗₂)/2.

"""
χₛ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2


"""
    χₐ(M₁, M₂, χ⃗₁, χ⃗₂)

Antisymmetric spin vector (χ⃗₁-χ⃗₂)/2.

"""
χₐ(M₁, M₂, χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
