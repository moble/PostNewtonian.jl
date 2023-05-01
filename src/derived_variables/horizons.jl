"""
    rₕ₁(s)

Horizon radius of black hole 1.

As defined on page 2, line 4, of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
rₕ₁(s::VecOrPNSystem) = M₁(s) * (1 + √(1-min(χ₁²(s),1)))

"""
    rₕ₂(s)

Horizon radius of black hole 2.

As defined on page 2, line 4, of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
rₕ₂(s::VecOrPNSystem) = M₂(s) * (1 + √(1-min(χ₂²(s),1)))


"""
    Ωₕ₁(s)

Horizon angular velocity of black hole 1.

As defined on page 2, line 5, of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
Ωₕ₁(s::VecOrPNSystem) = χ₁(s) / 2rₕ₁(s)

"""
    Ωₕ₂(s)

Horizon angular velocity of black hole 2.

As defined on page 2, line 5, of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
Ωₕ₂(s::VecOrPNSystem) = χ₂(s) / 2rₕ₂(s)


"""
    sin²θ₁(s)

Sine-squared of angle between spin of black hole 1 and vector to black hole 2.

Compare to Eq. (18) of [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function sin²θ₁(s::VecOrPNSystem)
    let χ₁²=χ₁²(s)
        ifelse(iszero(χ₁²), one(χ₁²), abs2vec(n̂(s)×χ⃗₁(s))/χ₁²)
    end
end

"""
    sin²θ₂(s)

Sine-squared of angle between spin of black hole 2 and vector to black hole 1.

Compare to Eq. (18) of [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function sin²θ₂(s::VecOrPNSystem)
    let χ₂²=χ₂²(s)
        ifelse(iszero(χ₂²), one(χ₂²), abs2vec(n̂(s)×χ⃗₂(s))/χ₂²)
    end
end


@doc raw"""
    ϕ̇̂₁(s)

Rate of rotation of black hole 2 about the spin of black hole 1, relative to orbital
rotation rate.

This is the rotation rate divided by ``v^3  = M\, \Omega``.

As defined in Eq. (19) of [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function ϕ̇̂₁(s::VecOrPNSystem)
    let χ₁=χ₁(s), sin²θ₁=sin²θ₁(s), M=M(s)
        ifelse(
            iszero(χ₁),
            inv(M),
            ifelse(
                iszero(sin²θ₁),
                zero(χ₁),
                ℓ̂(s) ⋅ χ⃗₁(s) / (M*χ₁*sin²θ₁)
            )
        )
    end
end

@doc raw"""
    ϕ̇̂₂(s)

Rate of rotation of black hole 1 about the spin of black hole 2, relative to orbital
rotation rate.

This is the rotation rate divided by ``v^3  = M\, \Omega``.

As defined in Eq. (19) of [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).
See the documentation section on ["Horizons"](@ref Horizons) for more details.
"""
function ϕ̇̂₂(s::VecOrPNSystem)
    let χ₂=χ₂(s), sin²θ₂=sin²θ₂(s), M=M(s)
        ifelse(
            iszero(χ₂),
            inv(M),
            ifelse(
                iszero(sin²θ₂),
                zero(χ₂),
                ℓ̂(s) ⋅ χ⃗₂(s) / (M*χ₂*sin²θ₂)
            )
        )
    end
end


@doc raw"""
    Î₀₁(s)

Horizon moment of inertia of black hole 1.

This is the moment divided by ``v^{12}``, as given by Eq. (10) of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
function Î₀₁(s::VecOrPNSystem)
    let χ₁²=χ₁²(s), sin²θ₁=sin²θ₁(s)
        (16rₕ₁(s) / 5M(s)^2) * M₁(s)^3 * sin²θ₁ * (
            1 - 3//4 * χ₁² + 15//4 * χ₁² * sin²θ₁
        )
    end
end

@doc raw"""
    Î₀₂(s)

Horizon moment of inertia of black hole 2.

This is the moment divided by ``v^{12}``, as given by Eq. (10) of [Alvi
(2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  See the documentation section
on ["Horizons"](@ref Horizons) for more details.
"""
function Î₀₂(s::VecOrPNSystem)
    let χ₂²=χ₂²(s), sin²θ₂=sin²θ₂(s)
        (16rₕ₂(s) / 5M(s)^2) * M₂(s)^3 * sin²θ₂ * (
            1 - 3//4 * χ₂² + 15//4 * χ₂² * sin²θ₂
        )
    end
end
