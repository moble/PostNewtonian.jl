"""
    μ(M₁, M₂)
    reduced_mass(M₁, M₂)

Compute the reduced mass ``(M₁ M₂)/(M₁+M₂)``.

"""
μ(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)
const reduced_mass = μ


"""
    ν(M₁, M₂)
    reduced_mass_ratio(M₁, M₂)

Compute the reduced mass ratio ``(M₁ M₂)/(M₁+M₂)^2``.

Note that the denominator is squared, unlike in the reduced mass [`μ`](@ref).

"""
ν(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)^2
const reduced_mass_ratio = ν
ν(;q) = q / (1+q)^2


"""
    δ(M₁, M₂)
    mass_difference_ratio(M₁, M₂)

Compute mass-difference ratio ``(M₁-M₂)/(M₁+M₂)``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that
``δ`` always be positive (or always negative), you are responsible for ensuring
that.

"""
δ(M₁, M₂) = (M₁ - M₂) / (M₁ + M₂)
const mass_difference_ratio = δ
δ(;q) = (q-1) / (q+1)


"""
    q(M₁, M₂)
    mass_ratio(M₁, M₂)

Compute mass ratio ``M₁/M₂``.

Note that we do not restrict to ``M₁ ≥ M₂`` or vice versa; if you prefer that
``q`` always be greater than or equal to 1 (or vice versa), you are responsible
for ensuring that.

"""
q(M₁, M₂) = M₁ / M₂
const mass_ratio = q


"""
    ℳ(M₁, M₂)
    chirp_mass(M₁, M₂)

Compute the chirp mass ℳ, which determines the leading-order orbital evolution
of a binary system due to energy loss by gravitational-wave emission.

The chirp mass is defined as
```math
  \\mathcal{M} = \\frac{(M₁ M₂)^{3/5}} {(M₁ + M₂)^{1/5}}.
```

"""
function ℳ(M₁, M₂)
    ((M₁ * M₂)^3 / (M₁ + M₂))^(1//5)
end
const chirp_mass = ℳ
