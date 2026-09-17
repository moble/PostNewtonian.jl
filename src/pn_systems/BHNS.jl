"""
    BHNS{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There is an additional field `Λ₂`
holding the (constant) tidal-coupling parameter of the neutron star.

Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`, and `Λ₂` always
refer to it; `M₁` and `χ⃗₁` always refer to the black hole.  See also [`NSNS`](@ref).
"""
struct BHNS{NT,ST<:DenseVector{NT},PNOrder} <: QuasisphericalSystem{NT,ST,PNOrder}
    state::ST
    Λ₂::NT
    BHNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function BHNS(; PNOrder=max_pn_order, kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(BHNS; PNOrder=PNOrder, kwargs...)
        return new{NT,ST,PNOrder}(state)
    end
    function BHNS(state; PNOrder=max_pn_order)
        @assert length(state) == length(symbols(BHNS))
        return new{eltype(state),typeof(state),prepare_pn_order(PNOrder)}(state)
    end
end

# The following are methods of functions defined in `state_variables.jl`,
# specialized for `BHNS` systems.
function pack_state(::Type{BHNS}; M₁, M₂, χ⃗₁, χ⃗₂, R=Rotor(1), v, Φ=0, Λ₂)
    [M₁; M₂; vec(QuatVec(χ⃗₁)); vec(QuatVec(χ⃗₂)); components(Rotor(R)); v; Φ; Λ₂]
end

state(pnsystem::BHNS) = pnsystem.state

function symbols(::Type{<:BHNS})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ, :Λ₂)
end
function ascii_symbols(::Type{<:BHNS})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi, :Lambda2)
end
