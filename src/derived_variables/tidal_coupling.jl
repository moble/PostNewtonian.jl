@doc raw"""
    Λ̃(pnsystem)
    Lambda_tilde(pnsystem)

Effective tidal deformability of the system.

Much as the chirp mass is a particularly measurable combination of the masses of the two
components of the binary, this quantity is a particularly measurable combination of the
tidal couplings ``\Lambda_1`` and ``\Lambda_2`` of the two components.

[Raithel et al. (2018)](https://doi.org/10.3847/2041-8213/aabcbf) define this as
```math
\tilde{\Lambda} = \frac{16}{13}
    \frac{(M_1 + 12 M_2) M_1^4 \Lambda_1 + (M_2 + 12 M_1) M_2^4 \Lambda_2} {M^5}.
```

See also [`Λ₁`](@ref) and [`Λ₂`](@ref).
"""
Λ̃(s::VecOrPNSystem) =
    16//13 * ((M₁(s) + 12M₂(s)) * M₁(s)^4 * Λ₁(s) + (M₂(s) + 12M₁(s)) * M₂(s)^4 * Λ₂(s)) /
    M(s)^5
const Lambda_tilde = Λ̃
