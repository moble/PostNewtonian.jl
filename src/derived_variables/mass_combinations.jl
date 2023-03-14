"""
    M(M₁, M₂)
    total_mass(M₁, M₂)

Compute the total mass ``M₁+M₂``.
"""
M(M₁, M₂) = M₁ + M₂
M(s::PNSystem) = M(M₁(s), M₂(s))
const total_mass = M


"""
    μ(M₁, M₂)
    reduced_mass(M₁, M₂)

Compute the reduced mass ``(M₁ M₂)/(M₁+M₂)``.
"""
μ(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)
μ(s::PNSystem) = μ(M₁(s), M₂(s))
const reduced_mass = μ


"""
    ν(M₁, M₂)
    reduced_mass_ratio(M₁, M₂)

Compute the reduced mass ratio ``(M₁ M₂)/(M₁+M₂)^2``.

Note that the denominator is squared, unlike in the reduced mass [`μ`](@ref).
"""
ν(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)^2
ν(s::PNSystem) = ν(M₁(s), M₂(s))
ν(;q) = q / (1+q)^2
const reduced_mass_ratio = ν


"""
    δ(M₁, M₂)
    mass_difference_ratio(M₁, M₂)

Compute mass-difference ratio ``(M₁-M₂)/(M₁+M₂)``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that ``δ`` always
be positive (or always negative), you are responsible for ensuring that.
"""
δ(M₁, M₂) = (M₁ - M₂) / (M₁ + M₂)
δ(s::PNSystem) = δ(M₁(s), M₂(s))
δ(;q) = (q-1) / (q+1)
const mass_difference_ratio = δ


"""
    q(M₁, M₂)
    mass_ratio(M₁, M₂)

Compute mass ratio ``M₁/M₂``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that ``q`` always
be greater than or equal to 1 (or vice versa), you are responsible for ensuring that.
"""
q(M₁, M₂) = M₁ / M₂
q(s::PNSystem) = q(M₁(s), M₂(s))
const mass_ratio = q


"""
    ℳ(M₁, M₂)
    chirp_mass(M₁, M₂)

Compute the chirp mass ℳ, which determines the leading-order orbital evolution of a binary
system due to energy loss by gravitational-wave emission.

The chirp mass is defined as
```math
  \\mathcal{M} = \\frac{(M_1 M_2)^{3/5}} {(M_1 + M_2)^{1/5}}.
```
"""
ℳ(M₁, M₂) = ((M₁ * M₂)^3 / (M₁ + M₂))^(1//5)
ℳ(s::PNSystem) = ℳ(M₁(s), M₂(s))
const chirp_mass = ℳ
