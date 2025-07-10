@doc raw"""
Expressions from [Blanchet2024](@cite)

This package mostly derives its notation from that paper.  The main exception is that
Blanchet uses the symbol ``m`` for mass and ``\mathup{M}`` for ADM mass, while we reserve
``m`` for the azimuthal number of spin-weighted spherical harmonics, and use ``M`` for mass
instead.  (ADM mass does not appear in this package as of this writing.)

In Sec. 3.2.4, Blanchet defines the total mass ``m=m₁+m₂``, the relative mass difference
``\Delta = (m₁ - m₂) / m``, the reduced mass ``μ = m₁ m₂ / m``, and the symmetric mass ratio
``ν ≡ \mu/m = m₁ m₂ / m²``.  He also poses ``X₁ = m₁ / m`` and ``X₂ = m₂ / m`` so that ``Δ =
X₁ - X₂`` and ``ν = X₁ X₂``.

In Eq. (369), he defines
```math
\gamma = \frac{GM}{rc^2}.
```

In Eq. (375), he defines the "frequency-related parameter"
```math
x = \left(\frac{GM\Omega}{c^3}\right)^{2/3},
```
Note that this ``\Omega`` is "the orbital frequency [...] of the circular orbit as measured
by a distant observer."  He does not define this completely, except to say that it is the
frequency associated with the asymptotic Killing vector field "in any natural coordinate
system which respects the helical symmetry" for an idealized circular system.

[Trestini2025](@cite) describes roughly the same frequency as the "waveform frequency" Ω,
and uses "orbital frequency" ω to refer to — presumably — the coordinate frequency.

"""
@pn_reference module Blanchet2024

import PostNewtonian: G, c, x, ν, γₑ

"""
    Eq483(pnsystem)

Quasi-circular non-spinning 4.5PN energy flux terms.

As given in Eq. (483) of [Blanchet2024](@cite), this expression includes "contributions from
tails, iterated tails, and tails-of-memory", but no contributions from spin.

"""
@pn_expression function Eq483(pnsystem)
    ℱ =
        (32c^5 / 5G * ν * x^5) * @pn_expansion(
            1 +
                (-1247/336 - 35/12 * ν)x +
                (4π)x^(3/2) +
                (-44711/9072 + 9271/504 * ν + 65/18 * ν^2)x^2 +
                (-8191/672 - 583/24 * ν)π * x^(5/2) +
                (
                    6643739519/69854400 + 16/3 * π^2 / 3 - 1712/105 * γₑ -
                    856/105 * ln(16x) + (-134543/7776 + 41/48 * π^2)ν - 94403/3024 * ν^2 -
                    775/324 * ν^3
                )x^3 +
                (-16285/504 + 214745/1728 * ν + 193385/3024 * ν^2)π * x^(7/2) +
                (
                    -323105549467/3178375200 + 232597/4410 * γₑ - 1369/126 * π^2 +
                    39931/294 * ln(2) - 47385/1568 * ln(3) +
                    232597/8820 * ln(x) +
                    (
                        -1452202403629/1466942400 + 41478/245 * γₑ - 267127/4608 * π^2 +
                        479062/2205 * ln(2) +
                        47385/392 * ln(3) +
                        20739/245 * ln(x)
                    )ν +
                    (1607125/6804 - 3157/384 * π^2)ν^2 +
                    6875/504 * ν^3 +
                    5/6 * ν^4
                )x^4 +
                (
                    265978667519/745113600 - 6848/105 * γₑ - 3424/105 * ln(16x) +
                    (2062241/22176 + 41/12 * π^2)ν - 133112905/290304 * ν^2 -
                    3719141/38016 * ν^3
                )π * x^(9/2)
        )
end

end  # module Blanchet2024
