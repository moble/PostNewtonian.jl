"""
    rₕ₁(s)

Horizon radius of black hole 1.

As defined on page 2, line 4, of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
rₕ₁(s::VecOrPNSystem) = M₁(s) * (1 + √(1 - min(χ₁²(s), 1)))

"""
    rₕ₂(s)

Horizon radius of black hole 2.

As defined on page 2, line 4, of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
rₕ₂(s::VecOrPNSystem) = M₂(s) * (1 + √(1 - min(χ₂²(s), 1)))

"""
    Ωₕ₁(s)

Horizon angular velocity of black hole 1.

As defined on page 2, line 5, of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
Ωₕ₁(s::VecOrPNSystem) = χ₁(s) / 2rₕ₁(s)

"""
    Ωₕ₂(s)

Horizon angular velocity of black hole 2.

As defined on page 2, line 5, of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
Ωₕ₂(s::VecOrPNSystem) = χ₂(s) / 2rₕ₂(s)

"""
    sin²θ₁(s)

Sine-squared of angle between spin of black hole 1 and vector to black hole 2.

Compare to Eq. (18) of [Alvi (2001)](https://doi.org/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function sin²θ₁(s::VecOrPNSystem)
    let χ₁² = χ₁²(s)
        ifelse(iszero(χ₁²), one(χ₁²), abs2vec(n̂(s) × χ⃗₁(s)) / χ₁²)
    end
end

"""
    sin²θ₂(s)

Sine-squared of angle between spin of black hole 2 and vector to black hole 1.

Compare to Eq. (18) of [Alvi (2001)](https://doi.org/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function sin²θ₂(s::VecOrPNSystem)
    let χ₂² = χ₂²(s)
        ifelse(iszero(χ₂²), one(χ₂²), abs2vec(n̂(s) × χ⃗₂(s)) / χ₂²)
    end
end

@doc raw"""
    ϕ̇̂₁(s)

Rate of rotation of black hole 2 about the spin of black hole 1, relative to orbital
rotation rate.

This is the rotation rate ϕ̇ as defined in Eq. (19) of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020), divided by ``v^3  = M\,
\Omega``.  This division is done to make sure we can track the relative PN order of terms
that depend on this.

See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function ϕ̇̂₁(s::VecOrPNSystem)
    let χ₁ = χ₁(s), sin²θ₁ = sin²θ₁(s), M = M(s)
        ifelse(
            iszero(χ₁),
            inv(M),
            ifelse(iszero(sin²θ₁), zero(χ₁), ℓ̂(s) ⋅ χ⃗₁(s) / (M * χ₁ * sin²θ₁)),
        )
    end
end

@doc raw"""
    ϕ̇̂₂(s)

Rate of rotation of black hole 1 about the spin of black hole 2, relative to orbital
rotation rate.

This is the rotation rate ϕ̇ as defined in Eq. (19) of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020), divided by ``v^3  = M\,
\Omega``.  This division is done to make sure we can track the relative PN order of terms
that depend on this.

See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function ϕ̇̂₂(s::VecOrPNSystem)
    let χ₂ = χ₂(s), sin²θ₂ = sin²θ₂(s), M = M(s)
        ifelse(
            iszero(χ₂),
            inv(M),
            ifelse(iszero(sin²θ₂), zero(χ₂), ℓ̂(s) ⋅ χ⃗₂(s) / (M * χ₂ * sin²θ₂)),
        )
    end
end

@doc raw"""
    Î₀₁(s)

Horizon moment of inertia of black hole 1.

This is the moment divided by ``ν^2 v^{12}``, as given by Eq. (10) of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
function Î₀₁(s::VecOrPNSystem)
    let χ₁² = χ₁²(s), sin²θ₁ = sin²θ₁(s)
        (16rₕ₁(s) / 5M(s)^2) * M₁(s)^3 * sin²θ₁ * (1 - 3//4 * χ₁² + 15//4 * χ₁² * sin²θ₁)
    end
end

@doc raw"""
    Î₀₂(s)

Horizon moment of inertia of black hole 2.

This is the moment divided by ``ν^2 v^{12}``, as given by Eq. (10) of [Alvi
(2001)](https://doi.org/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
function Î₀₂(s::VecOrPNSystem)
    let χ₂² = χ₂²(s), sin²θ₂ = sin²θ₂(s)
        (16rₕ₂(s) / 5M(s)^2) * M₂(s)^3 * sin²θ₂ * (1 - 3//4 * χ₂² + 15//4 * χ₂² * sin²θ₂)
    end
end

@doc raw"""
    κ₁(s)

The "quadrupolar polarisability" of object 1 used by [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).

Note that Bohé et al. refer to the closely related (and co-authored) [Marsat
(2014)](https://arxiv.org/abs/1411.4118), who notes above Eq. (4.7) that this is denoted
``C_{\mathrm{ES}^2}`` in [Levi and Steinhoff (2014)](https://arxiv.org/abs/1410.2601), who
in turn note that "we can set ... the Wilson coefficients ``C_{\mathrm{ES}^2} =
C_{\mathrm{BS}^3} = 1`` for the black hole case."

However, a very similar constant ``\kappa_A`` is used in Eq. (2.1) of [Buonanno et al.
(2012)](https://arxiv.org/abs/1209.6349).  They say that ``\kappa_A=1`` for an *isolated*
black hole, but can deviate from 1 for a black hole in a binary — though those deviations
occur at "much higher" order than those they consider (2PN).

See also [`λ₁`](@ref).

!!! warn
    This function will be incorrect for objects other than black holes.  It is not clear to
    me if this is the same quantity as ``C_Q`` used in some papers, such as [Bini and
    Geralico (2014)](https://arxiv.org/abs/1408.5261), but they point out that for neutron
    stars, the value varies between 4.3 and 7.4, depending on the equation of state.  This
    quantity may also be related to [`λ₁`](@ref).  Pull requests or issues with more
    information are welcome.
"""
function κ₁(s::VecOrPNSystem)
    return 1
end

"""
    κ₂(s)

The "quadrupolar polarisability" of object 2 used by [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).  See [`κ₁`](@ref) for more details.
"""
function κ₂(s::VecOrPNSystem)
    return 1
end

"""
    κ₊(s)

Equal to [`κ₁`](@ref)` + `[`κ₂`](@ref); defined below Eq. (3.28) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).
"""
function κ₊(s::VecOrPNSystem)
    return κ₁(s) + κ₂(s)
end

"""
    κ₋(s)

Equal to [`κ₁`](@ref)` - `[`κ₂`](@ref); defined below Eq. (3.28) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).
"""
function κ₋(s::VecOrPNSystem)
    return κ₁(s) - κ₂(s)
end

@doc raw"""
    λ₁(s)

The "quadrupolar polarisability" of object 1 used by [Marsat
(2014)](https://arxiv.org/abs/1411.4118), who notes above Eq. (4.11) that this is denoted
``C_{\mathrm{BS}^3}`` in [Levi and Steinhoff (2014)](https://arxiv.org/abs/1410.2601), who
in turn note that "we can set ... the Wilson coefficients ``C_{\mathrm{ES}^2} =
C_{\mathrm{BS}^3} = 1`` for the black hole case."

See also [`κ₁`](@ref).

!!! warn
    This function will be incorrect for objects other than black holes.  It is not clear to
    me if this is the same quantity as ``C_Q`` used in some papers, such as [Bini and
    Geralico (2014)](https://arxiv.org/abs/1408.5261), but they point out that for neutron
    stars, the value varies between 4.3 and 7.4, depending on the equation of state.  This
    quantity may also be related to [`λ₁`](@ref).  Pull requests or issues with more
    information are welcome.
"""
function λ₁(s::VecOrPNSystem)
    return 1
end

"""
    λ₂(s)

The "quadrupolar polarisability" of object 2 used by [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).  See [`λ₁`](@ref) for more details.
"""
function λ₂(s::VecOrPNSystem)
    return 1
end

"""
    λ₊(s)

Equal to [`λ₁`](@ref)` + `[`λ₂`](@ref); defined below Eq. (3.28) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).
"""
function λ₊(s::VecOrPNSystem)
    return λ₁(s) + λ₂(s)
end

"""
    λ₋(s)

Equal to [`λ₁`](@ref)` - `[`λ₂`](@ref); defined below Eq. (3.28) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).
"""
function λ₋(s::VecOrPNSystem)
    return λ₁(s) - λ₂(s)
end
