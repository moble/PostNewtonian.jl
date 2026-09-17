module FundamentalVariables

using ..PostNewtonian
using ..PostNewtonian: PNSystem, BBH, BHNS, NSNS, FDPNSystem, symbols
# using ..PostNewtonian: M₁index, M₂index, χ⃗₁indices, χ⃗₂indices, Rindices, vindex, Φindex
using Quaternionic: Quaternionic, QuatVec, Rotor

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₁, Λ₂, M1, M2, chi1, chi2, Phi, Lambda1, Lambda2, symbol_index

## NOTE:
## This indices used below are intimately bound to choices made in the definitions of
## the various `PNSystem`s.  Any changes there must be mirrored here, and vice versa.

"""
    M₁(pnsystem)
    M1(pnsystem)

Mass of object 1 in this system.
"""
function M₁(::T) where {T<:PNSystem}
    error("M₁ is not (yet) defined for PNSystem subtype `$T`.")
end
const M1 = M₁

"""
    M₂(pnsystem)
    M2(pnsystem)

Mass of object 2 in this system.
"""
function M₂(::T) where {T<:PNSystem}
    error("M₂ is not (yet) defined for PNSystem subtype `$T`.")
end
const M2 = M₂

"""
    χ⃗₁(pnsystem)
    chi1(pnsystem)

Dimensionless spin vector of object 1 in this system, as a `QuatVec`.
"""
χ⃗₁(s::T) where {T<:PNSystem} = @inbounds QuatVec(s[ :χ⃗₁ˣ : :χ⃗₁ᶻ ])

const chi1 = χ⃗₁

"""
    χ⃗₂(pnsystem)
    chi2(pnsystem)

Dimensionless spin vector of object 2 in this system, as a `QuatVec`.
"""
χ⃗₂(s::T) where {T<:PNSystem} = @inbounds QuatVec(s[ :χ⃗₂ˣ : :χ⃗₂ᶻ ])

const chi2 = χ⃗₂

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
R(s::T) where {T<:PNSystem} = @inbounds Rotor(s[ :Rʷ : :Rᶻ ])

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
function v(::T) where {T<:PNSystem}
    error("v is not (yet) defined for PNSystem subtype `$T`.")
end
v(; Ω, M=1) = ∛(M * Ω)

"""
    Φ(pnsystem)
    Phi(pnsystem)

Integrated orbital phase of the system.  It is computed as the integral of [`Ω`](@ref).
"""
function Φ(::T) where {T<:PNSystem}
    error("Φ is not (yet) defined for PNSystem subtype `$T`.")
end
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
quantity.  (`BHNS` systems can only have a nonzero value for ``\Lambda_2``.)  All other
types return `0`, which Julia can use to eliminate code that would then be 0.  Thus, it is
safe and efficient to use this quantity in any PN expression that specializes on the type of
`pnsystem`.

See also [`Λ₂`](@ref) and [`Λ̃`](@ref).
"""
Λ₁(pn::PNSystem) = zero(eltype(pn))
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

See also [`Λ₁`](@ref) and [`Λ̃`](@ref).
"""
Λ₂(pn::PNSystem) = zero(eltype(pn))
const Lambda2 = Λ₂

############################################################
# TODO This got moved here out of pn_systems.jl
for PNT ∈ (BBH, BHNS, NSNS)
    for (i, symbol) ∈ enumerate(symbols(PNT))
        # This will define, e.g., `M₁(pnsystem::BBH) = pnsystem.state[1]`.  We
        # could do this manually, but this is more concise and less error-prone.
        @eval PostNewtonian.FundamentalVariables begin
            $(symbol)(pnsystem::$PNT) = @inbounds pnsystem.state[$i]
            $(symbol)(pnsystem::FDPNSystem{NT,$PNT{NT,ST,PNOrder},PNOrder}) where {NT,ST,PNOrder} = @inbounds pnsystem.state[$i]
            function symbol_index(::Type{T}, ::Val{Symbol($symbol)}) where {T<:$PNT}
                $i
            end
            function symbol_index(::Type{FDPNSystem{NT,$PNT{NT,ST,PNOrder},PNOrder}}, ::Val{Symbol($symbol)}) where {NT,ST,PNOrder}
                $i
            end
        end
    end
end

symbol_index(t, s::Symbol) = symbol_index(t, Val(s))

end
