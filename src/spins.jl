"""
    χ⃗(S⃗, m)

Normalized spin vector ``S⃗/m²``.

"""
χ⃗(S⃗, m) = S⃗ / m^2


"""
    S(m₁, m₂, χ⃗₁, χ⃗₂)

Total spin vector ``S⃗₁+S⃗₂``.

"""
S(m₁, m₂, χ⃗₁, χ⃗₂) = χ⃗₁ * m₁^2 + χ⃗₂ * m₂^2


"""
    Σ(m₁, m₂, χ⃗₁, χ⃗₂)

Differential spin vector ``M(a⃗₂-a⃗₁)``.

"""
Σ(m₁, m₂, χ⃗₁, χ⃗₂) =  (m₁ + m₂) * (χ⃗₂ * m₂ - χ⃗₁ * m₁)


"""
    χₛ(m₁, m₂, χ⃗₁, χ⃗₂)

Symmetric spin vector (χ⃗₁+χ⃗₂)/2.

"""
χₛ(m₁, m₂, χ⃗₁, χ⃗₂) = (χ⃗₁ + χ⃗₂) / 2


"""
    χₐ(m₁, m₂, χ⃗₁, χ⃗₂)

Antisymmetric spin vector (χ⃗₁-χ⃗₂)/2.

"""
χₐ(m₁, m₂, χ⃗₁, χ⃗₂) = (χ⃗₁ - χ⃗₂) / 2
