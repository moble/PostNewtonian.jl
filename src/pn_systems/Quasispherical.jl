# NOTE: This type is documented separately, below, so that we can call
# `symbols(Quasispherical)` and `ascii_symbols(Quasispherical)` in the docstring.
@public abstract type Quasispherical{NT,PNOrder,ST} <: PNSystem{NT,PNOrder,ST} end

function symbols(::Type{<:Quasispherical})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end
function ascii_symbols(::Type{<:Quasispherical})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end

@doc """
    Quasispherical{NT,PNOrder,ST} <: PNSystem{NT,PNOrder,ST}

This is an abstract type for quasispherical post-Newtonian systems.  The built-in `BBH`,
`BHNS`, and `NSNS` systems subtype this.

The quantities common to all subtypes are

    $(symbols(Quasispherical))

which have the ASCII equivalents

    $(ascii_symbols(Quasispherical))

Note `BBH` has exactly these symbols, while `BHNS` and `NSNS` have additional symbols for
tidal-coupling parameters, `Λ₁` and `Λ₂` (or `Lambda1` and `Lambda2`).
"""
Quasispherical

function (::Type{T})(state, PNOrder=typemax(Int)) where {T<:Quasispherical}
    T{eltype(state),prepare_pn_order(PNOrder),typeof(state)}(state)
end

function prepare_Quasispherical(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=0, PNOrder=typemax(Int))
    χ⃗₁ = QuatVec(χ⃗₁)
    χ⃗₂ = QuatVec(χ⃗₂)
    state = MVector{14}(
        M₁, M₂, χ⃗₁[2], χ⃗₁[3], χ⃗₁[4], χ⃗₂[2], χ⃗₂[3], χ⃗₂[4], R[1], R[2], R[3], R[4], v, Φ
    )
    NT = eltype(state)
    PNOrder = prepare_pn_order(PNOrder)
    return (NT, PNOrder, state)
end
