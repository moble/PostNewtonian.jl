"""
    M(pnsystem)
    M(M₁, M₂)
    total_mass(pnsystem)
    total_mass(M1, M2)

Compute the total mass ``M₁+M₂``.
"""
M(M₁, M₂) = M₁ + M₂
M(s::VecOrPNSystem) = M(M₁(s), M₂(s))
const total_mass = M

"""
    μ(pnsystem)
    μ(M₁, M₂)
    reduced_mass(pnsystem)
    reduced_mass(M1, M2)

Compute the reduced mass ``(M₁ M₂)/(M₁+M₂)``.
"""
μ(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)
μ(s::VecOrPNSystem) = μ(M₁(s), M₂(s))
const reduced_mass = μ

"""
    ν(pnsystem)
    ν(M₁, M₂)
    reduced_mass_ratio(pnsystem)
    reduced_mass_ratio(M1, M2)

Compute the reduced mass ratio ``(M₁ M₂)/(M₁+M₂)^2``.

Note that the denominator is squared, unlike in the reduced mass [`μ`](@ref).
"""
ν(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)^2
ν(s::VecOrPNSystem) = ν(M₁(s), M₂(s))
ν(; q) = q / (1 + q)^2
const reduced_mass_ratio = ν

"""
    δ(pnsystem)
    δ(M₁, M₂)
    mass_difference_ratio(pnsystem)
    mass_difference_ratio(M1, M2)

Compute mass-difference ratio ``(M₁-M₂)/(M₁+M₂)``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that ``δ`` always
be positive (or always negative), you are responsible for ensuring that.
"""
δ(M₁, M₂) = (M₁ - M₂) / (M₁ + M₂)
δ(s::VecOrPNSystem) = δ(M₁(s), M₂(s))
δ(; q) = (q - 1) / (q + 1)
const mass_difference_ratio = δ

"""
    q(pnsystem)
    q(M₁, M₂)
    mass_ratio(pnsystem)
    mass_ratio(M1, M2)

Compute mass ratio ``M₁/M₂``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that ``q`` always
be greater than or equal to 1 (or vice versa), you are responsible for ensuring that.
"""
q(M₁, M₂) = M₁ / M₂
q(s::VecOrPNSystem) = q(M₁(s), M₂(s))
const mass_ratio = q

"""
    ℳ(pnsystem)
    ℳ(M₁, M₂)
    chirp_mass(pnsystem)
    chirp_mass(M1, M2)

Compute the chirp mass ℳ, which determines the leading-order orbital evolution of a binary
system due to energy loss by gravitational-wave emission.

The chirp mass is defined as
```math
  \\mathcal{M} = \\frac{(M_1 M_2)^{3/5}} {(M_1 + M_2)^{1/5}}.
```
"""
ℳ(M₁, M₂) = ((M₁ * M₂)^3 / (M₁ + M₂))^(1//5)
ℳ(s::VecOrPNSystem) = ℳ(M₁(s), M₂(s))
const chirp_mass = ℳ

"""
    X₁(pnsystem)
    X₁(M₁, M₂)
    X1(pnsystem)
    X1(M1, M2)

Compute the reduced *individual* mass ``M₁/(M₁+M₂)``.
"""
X₁(M₁, M₂) = M₁ / (M₁ + M₂)
X₁(s::VecOrPNSystem) = X₁(M₁(s), M₂(s))
const X1 = X₁

"""
    X₂(pnsystem)
    X₂(M₁, M₂)
    X2(pnsystem)
    X2(M1, M2)

Compute the reduced *individual* mass ``M₂/(M₁+M₂)``.
"""
X₂(M₁, M₂) = M₂ / (M₁ + M₂)
X₂(s::VecOrPNSystem) = X₂(M₁(s), M₂(s))
const X2 = X₂
