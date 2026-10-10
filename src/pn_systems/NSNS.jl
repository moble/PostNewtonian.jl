"""
    NSNS{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There are two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
struct NSNS{NT,ST<:DenseVector{NT},PNOrder} <: QuasisphericalSystem{NT,ST,PNOrder}
    state::ST

    NSNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function NSNS(; PNOrder=max_pn_order, kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(NSNS; PNOrder=PNOrder, kwargs...)
        return new{NT,ST,PNOrder}(state)
    end
    function NSNS(state; PNOrder=max_pn_order)
        @assert length(state) == length(symbols(NSNS))
        return new{eltype(state),typeof(state),prepare_pn_order(PNOrder)}(state)
    end
end

const BNS = NSNS

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `NSNS` systems.

function pack_state(::Type{NSNS}; M₁, M₂, χ⃗₁, χ⃗₂, R=Rotor(1), v, Φ=0, Λ₁, Λ₂)
    [M₁; M₂; vec(QuatVec(χ⃗₁)); vec(QuatVec(χ⃗₂)); components(Rotor(R)); v; Φ; Λ₁; Λ₂]
end

state(pnsystem::NSNS) = pnsystem.state
function symbols(::Type{<:NSNS})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ, :Λ₁, :Λ₂, )
end
function ascii_symbols(::Type{<:NSNS})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi, :Lambda1, :Lambda2,)
end
