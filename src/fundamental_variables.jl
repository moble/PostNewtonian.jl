module FundamentalVariables

using ..PostNewtonian
using ..PostNewtonian: PNSystem, BHNS, NSNS
using ..PostNewtonian: M₁index, M₂index, χ⃗₁indices, χ⃗₂indices, Rindices, vindex, Φindex
using Quaternionic

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, λ₁, λ₂,
    M1, M2, chi1, chi2, Phi, lambda1, lambda2

## NOTE:
## This indices used below are intimately bound to choices made in the definitions of
## the various `PNSystem`s.  Any changes there must be mirrored here, and vice versa.

"""
    M₁(pnsystem)
    M1(pnsystem)

Mass of object 1 in this system.
"""
M₁(s::PNSystem) = M₁(s.state)
M₁(state::AbstractVector) = @inbounds state[M₁index]
const M1 = M₁

"""
    M₂(pnsystem)
    M2(pnsystem)

Mass of object 2 in this system.
"""
M₂(s::PNSystem) = M₂(s.state)
M₂(state::AbstractVector) = @inbounds state[M₂index]
const M2 = M₂

"""
    χ⃗₁(pnsystem)
    chi1(pnsystem)

Dimensionless spin vector of object 1 in this system, as a `QuatVec`.
"""
χ⃗₁(s::PNSystem) = χ⃗₁(s.state)
χ⃗₁(state::AbstractVector) = @inbounds QuatVec(view(state, χ⃗₁indices)...)
const chi1 = χ⃗₁

"""
    χ⃗₂(pnsystem)
    chi2(pnsystem)

Dimensionless spin vector of object 2 in this system, as a `QuatVec`.
"""
χ⃗₂(s::PNSystem) = χ⃗₂(s.state)
χ⃗₂(state::AbstractVector) = @inbounds QuatVec(view(state, χ⃗₂indices)...)
const chi2 = χ⃗₂

"""
    R(pnsystem)

Orientation of the binary, as a `Rotor`.  This `Rotor` will rotate the `ẑ` vector to be
along the orbital angular velocity vector and the `x̂` vector to be along the separation
vector pointing from object 2 to object 1.

See also [`n̂`](@ref PostNewtonian.n̂), [`λ̂`](@ref PostNewtonian.λ̂), [`ℓ̂`](@ref
PostNewtonian.ℓ̂).
"""
R(s::PNSystem) = R(s.state)
R(state::AbstractVector) = @inbounds Rotor(view(state, Rindices)...)

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
v(s::PNSystem) = v(s.state)
v(state::AbstractVector) = @inbounds state[vindex]
v(;Ω, M=1) = ∛(M*Ω)

"""
    Φ(pnsystem)
    Phi(pnsystem)

Integrated orbital phase of the system.  It is computed as the integral of [`Ω`](@ref).
"""
Φ(s::PNSystem) = Φ(s.state)
Φ(state::AbstractVector) = @inbounds state[Φindex]
const Phi = Φ

@doc raw"""
    λ₁(pnsystem)
    lambda1(pnsystem)

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
const lambda1 = λ₁

@doc raw"""
    λ₂(pnsystem)
    lambda2(pnsystem)

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
const lambda2 = λ₂

end
