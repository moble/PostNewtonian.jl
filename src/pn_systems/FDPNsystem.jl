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
struct FDPNSystem{NT,PN<:PNSystem{NT},PNOrder} <: PNSystem{FDNode,Vector{FDNode},PNOrder}
    state::Vector{FDNode}
    constants::Vector{FDNode}

    function FDPNSystem(::Type{PN}, PNOrder=typemax(Int)) where {NT,PN<:PNSystem{NT}}
        return new{NT,prepare_pn_order(PNOrder)}(
            [FDNode(s) for s ∈ state_symbols(PN)], [FDNode(s) for s ∈ constants_symbols(PN)]
        )
    end
end

symbols(pnsystem::FDPNSystem{NT,PN}) where {NT,PN} = symbols(PN)

function symbol_index(pnsystem::FDPNSystem{NT,PN}, s::Symbol) where {NT,PN}
    symbol_index(PN, Val(s))
end

## TODO: See if this method is needed

## The old code had this, but I think it would probably just cause errors.  It might be
## relied upon in the functions where we take derivatives — 𝓔′code and γₚₙ₀′ — but even if
## so, maybe we could work around it with another function.
#Base.eltype(::FDPNSystem{FT}) where {FT} = FT

function characterize(pnsystem::FDPNSystem)
    return SVector{16,FastDifferentiation.Node}(
        pnsystem.state[1],
        pnsystem.state[2],
        pnsystem.state[3],
        pnsystem.state[4],
        pnsystem.state[5],
        pnsystem.state[6],
        pnsystem.state[7],
        pnsystem.state[8],
        pnsystem.state[9],
        pnsystem.state[10],
        pnsystem.state[11],
        pnsystem.state[12],
        pnsystem.state[13],
        pnsystem.state[14],
        Λ₁(pnsystem),
        Λ₂(pnsystem),
    )
end
