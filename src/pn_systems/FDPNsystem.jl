"""
    FDPNSystem{NT, PN, PNOrder} <: PNSystem{FDNode, Vector{FDNode}, PNOrder}

A `PNSystem` that contains information as variables from
[`FastDifferentiation.jl`](https://docs.juliahub.com/General/FastDifferentiation/stable/).

Note that this type also involves the type parameter `PN`, which is actually the type of a
`PNSystem`, and its type parameter `NT`, which will be the number type of actual numbers
that eventually get fed into (and will be passed out from) functions that use this system.

One important example of what this type is used for is computing the derivative of the
orbital binding energy, `𝓔′` — and in particular, for generating the corresponding function
method to apply to a given `PNSystem`.
"""
@export struct FDPNSystem{NT,PNOrder,PN<:Type{<:PNSystem{NT,PNOrder}}} <:
               PNSystem{FDNode,PNOrder,Vector{FDNode}}
    state::Vector{FDNode}

    FDPNSystem(pnsystem::PNSystem) = FDPNSystem(typeof(pnsystem))
    function FDPNSystem(::Type{PN}) where {NT,PNOrder,PN<:PNSystem{NT,PNOrder}}
        return new{NT,PNOrder,PN}([FDNode(s) for s ∈ symbols(PN)])
    end
end

symbols(pnsystem::FDPNSystem{NT,PNOrder,PN}) where {NT,PNOrder,PN} = symbols(PN)

function symbol_index(pnsystem::FDPNSystem{NT,PNOrder,PN}, s::Symbol) where {NT,PNOrder,PN}
    symbol_index(PN, Val(s))
end

## TODO: See if this method is needed

## The old code had this, but I think it would probably just cause errors.  It might be
## relied upon in the functions where we take derivatives — 𝓔′code and γₚₙ₀′ — but even if
## so, maybe we could work around it with another function.
#Base.eltype(::FDPNSystem{FT}) where {FT} = FT
