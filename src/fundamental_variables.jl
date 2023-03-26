module FundamentalVariables

using ..PostNewtonian
using ..PostNewtonian: PNSystem, BHNS, NSNS
using Quaternionic

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, λ₁, λ₂

"""
    M₁(pnsystem)

Mass of object 1 in this system.
"""
M₁(s::PNSystem) = @inbounds s.state[1]

"""
    M₂(pnsystem)

Mass of object 2 in this system.
"""
M₂(s::PNSystem) = @inbounds s.state[2]

"""
    χ⃗₁(pnsystem)

Dimensionless spin vector of object 1 in this system, as a `QuatVec`.
"""
χ⃗₁(s::PNSystem) = @inbounds QuatVec(view(s.state, 3:5)...)

"""
    χ⃗₂(pnsystem)

Dimensionless spin vector of object 2 in this system, as a `QuatVec`.
"""
χ⃗₂(s::PNSystem) = @inbounds QuatVec(view(s.state, 6:8)...)

"""
    R(pnsystem)

Orientation of the binary, as a `Rotor`.  This `Rotor` will rotate the `ẑ` vector to be
along the orbital angular velocity vector and the `x̂` vector to be along the separation
vector pointing from object 2 to object 1.

See also [`n̂`](@ref PostNewtonian.n̂), [`λ̂`](@ref PostNewtonian.λ̂), [`ℓ̂`](@ref
PostNewtonian.ℓ̂).
"""
R(s::PNSystem) = @inbounds Rotor(view(s.state, 9:12)...)

@doc raw"""
    v(pnsystem)
    v(;Ω, M=1)

Post-Newtonian velocity parameter.  This is related to the orbital angular frequency
``\Omega`` as
```math
v \coloneq (M\,\Omega)^{1/3},
```
where ``M`` is the total mass of the binary.

Note that if you want to pass the value ``Ω`` (rather than a `PNSystem`), you must pass it
as a keyword argument — as in `v(Ω=0.1)`.

See also [`Ω`](@ref).
"""
v(s::PNSystem) = @inbounds s.state[13]
v(;Ω, M=1) = ∛(M*Ω)

"""
    Φ(pnsystem)

Integrated orbital phase of the system.  Note that by default, this is not available;
additional arguments must be passed to ensure that `pnsystem` tracks this variable.  When it
is available, it is given by the integral of `Ω`.

Note that if it is not available, a `BoundsError` will be raised by this function.
"""
Φ(s::PNSystem) = s.state[14]  # NO @inbounds

@doc raw"""
    λ₁(pnsystem)

Tidal coupling parameter of object 1 in this system.

This is related to the Love number ``k_2`` and the star's radius ``R`` as
```math
\lambda_1 = \frac{2}{3} k_2 R^5.
```
There is another common parameter denoted
```math
\Lambda_1 = \frac{3}{2} \frac{M_1}{M_2} \lambda_1.
```

Note that — as of this writing — only `NSNS` systems can have a nonzero value for this
quantity.  All other types return `0`, which Julia can use to eliminate code that would then
be 0.  Thus, it is safe and efficient to use this quantity in any PN expression that
specializes on the type of `pnsystem`.
"""
λ₁(::PNSystem) = 0
λ₁(pn::NSNS) = pn.λ₁
λ₁(pn::SymbolicPNSystem) = pn.λ₁

@doc raw"""
    λ₂(pnsystem)

Tidal coupling parameter of object 2 in this system.

This is related to the Love number ``k_2`` and the star's radius ``R`` as
```math
\lambda_2 = \frac{2}{3} k_2 R^5.
```
There is another common parameter denoted
```math
\Lambda_2 = \frac{3}{2} \frac{M_2}{M_1} \lambda_2.
```

Note that — as of this writing — only `BHNS` and `NSNS` systems can have a nonzero value for
this quantity.  All other types return `0`, which Julia can use to eliminate code that would
then be 0.  Thus, it is safe and efficient to use this quantity in any PN expression that
specializes on the type of `pnsystem`.
"""
λ₂(::PNSystem) = 0
λ₂(pn::BHNS) = pn.λ₂
λ₂(pn::NSNS) = pn.λ₂
λ₂(pn::SymbolicPNSystem) = pn.λ₂

end
