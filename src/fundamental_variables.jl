module FundamentalVariables

using ..PostNewtonian
using ..PostNewtonian: PNSystem, BHNS, NSNS
using ..PostNewtonian: M₁index, M₂index, χ⃗₁indices, χ⃗₂indices, Rindices, vindex, Φindex
using Quaternionic

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₁, Λ₂,
    M1, M2, chi1, chi2, Phi, Lambda1, Lambda2

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
v \colonequals (M\,\Omega)^{1/3},
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
quantity.  All other types return `0`, which Julia can use to eliminate code that would then
be 0.  Thus, it is safe and efficient to use this quantity in any PN expression that
specializes on the type of `pnsystem`.
"""
Λ₁(::PNSystem) = 0
Λ₁(pn::NSNS) = pn.Λ₁
Λ₁(pn::SymbolicPNSystem) = pn.Λ₁
const Lambda1 = Λ₁

@doc raw"""
    Λ₂(pnsystem)
    Lambda2(pnsystem)

Quadrupolar tidal coupling parameter of object 2 in this system.

See [`Λ₁`](@ref) for details about the definition, swapping "object 1" with "object 2".

Note that — as of this writing — only `BHNS` and `NSNS` systems can have a nonzero value for
this quantity.  All other types return `0`, which Julia can use to eliminate code that would
then be 0.  Thus, it is safe and efficient to use this quantity in any PN expression that
specializes on the type of `pnsystem`.
"""
Λ₂(::PNSystem) = 0
Λ₂(pn::BHNS) = pn.Λ₂
Λ₂(pn::NSNS) = pn.Λ₂
Λ₂(pn::SymbolicPNSystem) = pn.Λ₂
const Lambda2 = Λ₂

end
