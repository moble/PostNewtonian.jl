"""
    G(pnsystem)

Return Newton's gravitational constant for the given `pnsystem`.

By default, the value is one *with the same number type as `pnsystem`*.  It can be
overridden for subtypes of `PNSystem` that use different units or conventions.

However, note that this function should specialize on the number type of `pnsystem`, rather
than just returning the integer `1`, because there may be expressions with factors such as
`G/3` which will immediately convert to `Float64` if `G` is just `1`, so the result will not
have the expected precision.
"""
@public G(::PNSystem{NT}) where {NT} = one(NT)
G(::FDPNSystem{NT,PN}) where {NT,PN} = one(PN)

"""
    c(pnsystem)

Return the speed of light for the given `pnsystem`.

By default, the value is one *with the same number type as `pnsystem`*.  It can be
overridden for subtypes of `PNSystem` that use different units or conventions.

However, note that this function should specialize on the number type of `pnsystem`, rather
than just returning the integer `1`, because there may be expressions with factors such as
`c/3` which will immediately convert to `Float64` if `c` is just `1`, so the result will not
have the expected precision.
"""
@public c(::PNSystem{NT}) where {NT} = one(NT)
c(::FDPNSystem{NT,PN}) where {NT,PN} = one(PN)

"""
    M₁(pnsystem)
    M1(pnsystem)

Mass of object 1 in this system.
"""
@public function M₁(::T) where {T<:PNSystem}
    error("M₁ is not (yet) defined for PNSystem subtype `$T`.")
end
M₁(fdpnsystem::FDPNSystem) = fdpnsystem[:M₁]
@public const M1 = M₁

"""
    M₂(pnsystem)
    M2(pnsystem)

Mass of object 2 in this system.
"""
@public function M₂(::T) where {T<:PNSystem}
    error("M₂ is not (yet) defined for PNSystem subtype `$T`.")
end
M₂(fdpnsystem::FDPNSystem) = fdpnsystem[:M₂]
@public const M2 = M₂

"""
    χ₁ˣ(pnsystem)
    chi1x(pnsystem)

`x`-component of dimensionless spin vector of object 1 in this system, as a `QuatVec`.

See [`χ⃗₁`](@ref) for details.
"""
@public function χ₁ˣ(::T) where {T<:PNSystem}
    error("χ₁ˣ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₁ˣ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₁ˣ]
@public const chi1x = χ₁ˣ

"""
    χ₁ʸ(pnsystem)
    chi1y(pnsystem)

`y`-component of dimensionless spin vector of object 1 in this system, as a `QuatVec`.

See [`χ⃗₁`](@ref) for details.
"""
@public function χ₁ʸ(::T) where {T<:PNSystem}
    error("χ₁ʸ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₁ʸ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₁ʸ]
@public const chi1y = χ₁ʸ

"""
    χ₁ᶻ(pnsystem)
    chi1z(pnsystem)

`z`-component of dimensionless spin vector of object 1 in this system, as a `QuatVec`.

See [`χ⃗₁`](@ref) for details.
"""
@public function χ₁ᶻ(::T) where {T<:PNSystem}
    error("χ₁ᶻ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₁ᶻ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₁ᶻ]
@public const chi1z = χ₁ᶻ

"""
    χ₂ˣ(pnsystem)
    chi2x(pnsystem)

`x`-component of dimensionless spin vector of object 2 in this system, as a `QuatVec`.

See [`χ⃗₂`](@ref) for details.
"""
@public function χ₂ˣ(::T) where {T<:PNSystem}
    error("χ₂ˣ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₂ˣ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₂ˣ]
@public const chi2x = χ₂ˣ

"""
    χ₂ʸ(pnsystem)
    chi2y(pnsystem)

`y`-component of dimensionless spin vector of object 2 in this system, as a `QuatVec`.

See [`χ⃗₂`](@ref) for details.
"""
@public function χ₂ʸ(::T) where {T<:PNSystem}
    error("χ₂ʸ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₂ʸ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₂ʸ]
@public const chi2y = χ₂ʸ

"""
    χ₂ᶻ(pnsystem)
    chi2z(pnsystem)

`z`-component of dimensionless spin vector of object 2 in this system, as a `QuatVec`.

See [`χ⃗₂`](@ref) for details.
"""
@public function χ₂ᶻ(::T) where {T<:PNSystem}
    error("χ₂ᶻ is not (yet) defined for PNSystem subtype `$T`.")
end
χ₂ᶻ(fdpnsystem::FDPNSystem) = fdpnsystem[:χ₂ᶻ]
@public const chi2z = χ₂ᶻ

"""
    Rʷ(pnsystem)
    Rw(pnsystem)

Scalar component of the orientation `Rotor` of the binary.

See [`R`](@ref) for details.
"""
@public function Rʷ(::T) where {T<:PNSystem}
    error("Rʷ is not (yet) defined for PNSystem subtype `$T`.")
end
Rʷ(fdpnsystem::FDPNSystem) = fdpnsystem[:Rʷ]
@public const Rw = Rʷ

"""
    Rˣ(pnsystem)
    Rx(pnsystem)

`x`-component of the orientation `Rotor` of the binary.

See [`R`](@ref) for details.
"""
@public function Rˣ(::T) where {T<:PNSystem}
    error("Rˣ is not (yet) defined for PNSystem subtype `$T`.")
end
Rˣ(fdpnsystem::FDPNSystem) = fdpnsystem[:Rˣ]
@public const Rx = Rˣ

"""
    Rʸ(pnsystem)
    Ry(pnsystem)

`y`-component of the orientation `Rotor` of the binary.

See [`R`](@ref) for details.
"""
@public function Rʸ(::T) where {T<:PNSystem}
    error("Rʸ is not (yet) defined for PNSystem subtype `$T`.")
end
Rʸ(fdpnsystem::FDPNSystem) = fdpnsystem[:Rʸ]
@public const Ry = Rʸ

"""
    Rᶻ(pnsystem)
    Rz(pnsystem)

`z`-component of the orientation `Rotor` of the binary.

See [`R`](@ref) for details.
"""
@public function Rᶻ(::T) where {T<:PNSystem}
    error("Rᶻ is not (yet) defined for PNSystem subtype `$T`.")
end
Rᶻ(fdpnsystem::FDPNSystem) = fdpnsystem[:Rᶻ]
@public const Rz = Rᶻ

@doc raw"""
    v(pnsystem)
    v(;Ω, M=1)

Post-Newtonian velocity parameter.  This is related to the orbital angular frequency
``\Omega`` as
```math
v \colonequals (M\,\Omega)^{1/3},
```
where ``M`` is the total mass of the binary.

Note that if you want to pass the value ``Ω`` (rather than a `PNSystem`), you must pass it
as a keyword argument — as in `v(Ω=0.1)`.

See also [`Ω`](@ref).
"""
@public function v(::T) where {T<:PNSystem}
    error("v is not (yet) defined for PNSystem subtype `$T`.")
end
v(fdpnsystem::FDPNSystem) = fdpnsystem[:v]
v(; Ω, M=1) = ∛(M * Ω)

"""
    Φ(pnsystem)
    Phi(pnsystem)

Integrated orbital phase of the system.  It is computed as the integral of [`Ω`](@ref).
"""
@public function Φ(::T) where {T<:PNSystem}
    error("Φ is not (yet) defined for PNSystem subtype `$T`.")
end
Φ(fdpnsystem::FDPNSystem) = fdpnsystem[:Φ]
@public const Phi = Φ

@doc raw"""
    Λ₁(pnsystem)
    Lambda1(pnsystem)

Quadrupolar tidal-coupling parameter of object 1 in this system.

We imagine object 1 begin placed in an (adiabatic) external field with Newtonian potential
``\phi``, resulting in a tidal field measured by ``\partial_i \partial_j \phi`` evaluated at
the center of mass of the object.  This induces a quadrupole moment ``Q_{ij}`` in object 1,
which can be related to the tidal field as
```math
Q_{ij} = -\frac{G^4}{c^{10}} \Lambda_1 M_1^5 \partial_i \partial_j \phi,
```
where ``M_1`` is the mass of object 1.  This tidal-coupling parameter ``\Lambda_1`` can be
related to the Love number ``k_2`` (where the subscript 2 refers to the fact that this is
for the ``\ell=2`` quadrupole, rather than object 2) as
```math
\Lambda_1 = \frac{2}{3} \frac{c^{10}}{G^5} \frac{R_1^5}{M_1^5} k_2,
```
where ``R_1`` is the radius of object 1.  Note that ``\Lambda_1`` is dimensionless.  For
black holes, it is precisely zero; for neutron stars it may range up to 1; more exotic
objects may have significantly larger values.

Note that — as of this writing — only `NSNS` systems can have a nonzero value for this
quantity.  (`BHNS` systems can only have a nonzero value for ``\Lambda_2``.)  All other
types return `0`, which Julia can use to eliminate code that would then be 0.  Thus, it is
safe and efficient to use this quantity in any PN expression that specializes on the type of
`pnsystem`.

See also [`Λ₂`](@ref) and [`Λ̃`](@ref).
"""
@public function Λ₁(::T) where {T<:PNSystem}
    error("Λ₁ is not (yet) defined for PNSystem subtype `$T`.")
end
Λ₁(fdpnsystem::FDPNSystem) = fdpnsystem[:Λ₁]
@public const Lambda1 = Λ₁

@doc raw"""
    Λ₂(pnsystem)
    Lambda2(pnsystem)

Quadrupolar tidal coupling parameter of object 2 in this system.

See [`Λ₁`](@ref) for details about the definition, swapping "object 1" with "object 2".

Note that — as of this writing — only `BHNS` and `NSNS` systems can have a nonzero value for
this quantity.  All other types return `0`, which Julia can use to eliminate code that would
then be 0.  Thus, it is safe and efficient to use this quantity in any PN expression that
specializes on the type of `pnsystem`.

See also [`Λ₁`](@ref) and [`Λ̃`](@ref).
"""
@public function Λ₂(::T) where {T<:PNSystem}
    error("Λ₂ is not (yet) defined for PNSystem subtype `$T`.")
end
Λ₂(fdpnsystem::FDPNSystem) = fdpnsystem[:Λ₂]
@public const Lambda2 = Λ₂

#################################################################
# Not actually state variables, but aggregates of state variables

"""
    χ⃗₁(pnsystem)
    chi1(pnsystem)

Dimensionless spin vector of object 1 in this system, as a `QuatVec`.

See also [`χ₁ˣ`](@ref), [`χ₁ʸ`](@ref), and [`χ₁ᶻ`](@ref) for the individual components.
"""
@public function χ⃗₁(::T) where {T<:PNSystem}
    QuatVec(χ₁ˣ(pnsystem), χ₁ʸ(pnsystem), χ₁ᶻ(pnsystem))
end
@public const chi1 = χ⃗₁

"""
    χ⃗₂(pnsystem)
    chi2(pnsystem)

Dimensionless spin vector of object 2 in this system, as a `QuatVec`.

See also [`χ₂ˣ`](@ref), [`χ₂ʸ`](@ref), and [`χ₂ᶻ`](@ref) for the individual components.
"""
@public function χ⃗₂(::T) where {T<:PNSystem}
    QuatVec(χ₂ˣ(pnsystem), χ₂ʸ(pnsystem), χ₂ᶻ(pnsystem))
end
@public const chi2 = χ⃗₂

"""
    R(pnsystem)

Orientation of the binary, as a `Rotor`.

At any instant, the binary is represented by the right-handed triad ``(n̂, λ̂, ℓ̂)``, where
[``n̂``](@ref PostNewtonian.n̂) is the unit vector pointing from object 2 to object 1, and
the instantaneous velocities of the binary's elements are in the ``n̂``-``λ̂`` plane.  This
`Rotor` will rotate the ``x̂`` vector to be along ``n̂``,  the ``ŷ`` vector to be along
``λ̂``, and  the ``ẑ`` vector to be along ``ℓ̂``.

Note that the angular velocity associated to `R` is given by ``Ω⃗ = 2 Ṙ R̄ = Ω ℓ̂ + ϖ n̂``.
(Any component of ``Ω⃗`` along ``λ̂`` would violate the condition that the velocities be in
the ``n̂``-``λ̂`` plane.)  Here, the scalar quantity ``Ω`` is the orbital angular frequency,
and ``ϖ`` is the precession angular frequency.

See also [`n̂`](@ref PostNewtonian.n̂), [`λ̂`](@ref PostNewtonian.λ̂), [`ℓ̂`](@ref
PostNewtonian.ℓ̂), [`Ω`](@ref PostNewtonian.Ω), and [`𝛡`](@ref PostNewtonian.𝛡)``=ϖ n̂``.
"""
@public function R(pnsystem::T) where {NT,T<:PNSystem{NT}}
    # We use this explicit constructor (with type parameter) to avoid normalization
    # that would probably just complicate derivatives.
    Rotor{NT}(Rʷ(pnsystem), Rˣ(pnsystem), Rʸ(pnsystem), Rᶻ(pnsystem))
end
R(fdpnsystem::FDPNSystem) = fdpnsystem[:R]

@testitem "State variables" begin
    # Test BBH, BHNS, and NSNS state variables
    for pnsystem ∈ (BBH(randn(14)), BHNS(randn(15)), NSNS(randn(16)))
        for (i, (s, a)) ∈ enumerate(zip(symbols(pnsystem), ascii_symbols(pnsystem)))
            @test PostNewtonian.eval(s)(pnsystem) == pnsystem.state[i]
            @test PostNewtonian.eval(a)(pnsystem) == pnsystem.state[i]
        end
    end
    bbh = BBH(randn(14))
    @test PostNewtonian.Λ₁(bbh) == 0
    @test PostNewtonian.Λ₂(bbh) == 0
    @test PostNewtonian.Lambda1(bbh) == 0
    @test PostNewtonian.Lambda2(bbh) == 0
    bhns = BHNS(randn(15))
    @test PostNewtonian.Λ₁(bhns) == 0
    @test PostNewtonian.Lambda1(bhns) == 0

    # Test FDPNSystem state variables
    for pnsystem ∈ (BBH(randn(14)), BHNS(randn(15)), NSNS(randn(16)))
        fdpnsystem = FDPNSystem(pnsystem)
        for (i, (s, a)) ∈ enumerate(zip(symbols(pnsystem), ascii_symbols(pnsystem)))
            @test PostNewtonian.eval(s)(fdpnsystem).node_value == s
            @test PostNewtonian.eval(a)(fdpnsystem).node_value == s
        end
    end
end
