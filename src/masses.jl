"""
    μ(M₁, M₂)
    reduced_mass(M₁, M₂)

Compute the reduced mass ``(M₁ M₂)/(M₁+M₂)``.

"""
μ(M₁, M₂) = reduced_mass(M₁, M₂)
reduced_mass(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)


"""
    ν(M₁, M₂)
    reduced_mass_ratio(M₁, M₂)

Compute the reduced mass ratio ``(M₁ M₂)/(M₁+M₂)^2``.

Note that the denominator is squared, unlike in the reduced mass [`μ`](@ref).

"""
ν(M₁, M₂) = reduced_mass_ratio(M₁, M₂)
reduced_mass_ratio(M₁, M₂) = (M₁ * M₂) / (M₁ + M₂)^2
ν(q) = q / (1+q)^2


"""
    δ(M₁, M₂)
    mass_difference_ratio(M₁, M₂)

Compute mass-difference ratio ``(M₁-M₂)/(M₁+M₂)``.

"""
δ(M₁, M₂) = reduced_mass(M₁, M₂)
mass_difference_ratio(M₁, M₂) = (M₁ - M₂) / (M₁ + M₂)
δ(q) = (q-1) / (q+1)


"""
    q(M₁, M₂)
    mass_ratio(M₁, M₂)

Compute mass ratio ``M₁/M₂``.

Note that we do not restrict to `M₁ ≥ M₂` or vice versa; if you prefer that
``q`` always be greater than or equal to 1 (or vice versa), you are responsible
for ensuring that

"""
q(M₁, M₂) = mass_ratio(M₁, M₂)
mass_ratio(M₁, M₂) = M₁ / M₂


"""
    ℳ(M₁, M₂)
    chirp_mass(M₁, M₂)

Compute the chirp mass ℳ, which determines the leading-order orbital evolution
of a binary system due to energy loss by gravitational-wave emission.

The chirp mass is defined as
```math
  \\mathcal{M} = \frac{(M₁ M₂)^{3/5}} {(M₁ + M₂)^{1/5}}.
```

"""
function chirp_mass(M₁, M₂)
    T = promote_type(typeof(M₁), typeof(M₂))
    ((M₁ * M₂)^3 / (M₁ + M₂))^inv(T(5))
end
ℳ(M₁, M₂) = chirp_mass(M₁, M₂)
