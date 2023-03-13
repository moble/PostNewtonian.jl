abstract type AbstractPNSystem{T, PNOrder, Expansion} end

eltype(::AbstractPNSystem{T}) where {T} = T
pn_order(::AbstractPNSystem{T, PNOrder})::Irrational{Int} where {T, PNOrder} = PNOrder
expansion_type(::AbstractPNSystem{T, P, Expansion}) where {T, P, Expansion} = Expansion

struct BBH{T, PNOrder, Expansion} <: AbstractPNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
end
function BBH(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1
)
    state = [M₁; M₂; vec(χ⃗₁); vec(χ⃗₂); components(R); v]
    if !isnothing(Φ)
        state = [state; Φ]
    end
    T = eltype(state)
    PNOrder = if PNOrder!=typemax(Int)
        round(Int, 2PNOrder) // 2
    else
        typemax(Int) // 1
    end
    BBH{T, PNOrder, Expansion}(state)
end
