"""
    FDPNSystem{NT, PNOrder}(state, Λ₁, Λ₂)

TODO UPDATE
A `PNSystem` that contains information as variables from
[`FastDifferentiation.jl`](https://docs.juliahub.com/General/FastDifferentiation/stable/).

See also [`fd_pnsystem`](@ref) for a particular instance of this type.  Note that this type
also involves the type `NT`, which will be the numeric type of actual numbers that eventually
get fed into (and will be passed out from) functions that use this system.  The correct type
of `FDPNSystem` is used in calculating `𝓔′`.
"""
struct FDPNSystem{NT,PN<:PNSystem{NT},PNOrder} <: PNSystem{FastDifferentiation.Node,Vector{FastDifferentiation.Node},PNOrder}
    state::Vector{FastDifferentiation.Node}

    function FDPNSystem(::Type{PN}, PNOrder=max_pn_order) where {NT,PN<:PNSystem{NT}}
        return new{NT,PN,prepare_pn_order(PNOrder)}([FastDifferentiation.Node(s) for s ∈ symbols(PN)])
    end
end

state(pnsystem::FDPNSystem) = pnsystem.state

"""
    fd_pnsystem

A symbolic `PNSystem` that contains symbolic information for all types of `PNSystem`s.

In particular, note that this object has (essentially) infinite `PNOrder`, has nonzero
values for quantities like `Λ₁` and `Λ₂`, and assumes that the eventual output will be in
`Float64`.  If you want different choices, you may need to call [`FDPNSystem`](@ref)
yourself, or even construct a different specialized subtype of `PNSystem` (it's not hard).

# Examples
```jldoctest
julia> using PostNewtonian: M₁, M₂, χ⃗₁, χ⃗₂, FDPNSystem

julia> fd_pnsystem = FDPNSystem(Float64)
FDPNSystem{Float64, 9223372036854775805//2}(FastDifferentiation.Node[M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ], Λ₁, Λ₂)

julia> M₁(fd_pnsystem), M₂(fd_pnsystem)
(M₁, M₂)

julia> χ⃗₁(fd_pnsystem)
 + χ⃗₁ˣ𝐢 + χ⃗₁ʸ𝐣 + χ⃗₁ᶻ𝐤

julia> χ⃗₂(fd_pnsystem)
 + χ⃗₂ˣ𝐢 + χ⃗₂ʸ𝐣 + χ⃗₂ᶻ𝐤
```
"""
#const fd_pnsystem = FDPNSystem(Float64)

function StaticArrays.SVector(pnsystem::FDPNSystem)
    return SVector{length(pnsystem), FastDifferentiation.Node}(pnsystem.state)
end
