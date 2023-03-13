abstract type AbstractPNSystem{T, PNOrder, Expansion} end

eltype(::AbstractPNSystem{T}) where {T} = T
pn_order(::AbstractPNSystem{T, PNOrder})::Irrational{Int} where {T, PNOrder} = PNOrder
expansion_type(::AbstractPNSystem{T, P, Expansion}) where {T, P, Expansion} = Expansion

order_index(pn::AbstractPNSystem) = 1 + Int(2pn_order(pn))

function prepare_system(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1, kwargs...
)
    state = [M₁; M₂; vec(χ⃗₁); vec(χ⃗₂); components(R); v]
    if !isnothing(Φ)
        state = [state; Φ]
    end
    T = eltype(state)
    PNOrder = if PNOrder!=typemax(Int)
        round(Int, 2PNOrder) // 2
    else
        typemax(Int) // 2
    end
    (T, PNOrder, Expansion, state, kwargs)
end


struct BBH{T, PNOrder, Expansion} <: AbstractPNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
end
function BBH(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1
)
    (T, PNOrder, Expansion, state, kwargs) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    BBH{T, PNOrder, Expansion}(state)
end


struct BHNS{T, PNOrder, Expansion} <: AbstractPNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
    λ₂::T
end
function BHNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, λ₂, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1
)
    (T, PNOrder, Expansion, state, kwargs) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    BHNS{T, PNOrder, Expansion}(state, λ₂)
end
λ₂(pn::BHNS) = pn.λ₂


struct NSNS{T, PNOrder, Expansion} <: AbstractPNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
    λ₁::T
    λ₂::T
end
function NSNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, λ₁, λ₂, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1
)
    (T, PNOrder, Expansion, state, kwargs) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    NSNS{T, PNOrder, Expansion}(state, λ₁, λ₂)
end
λ₁(pn::NSNS) = pn.λ₁
λ₂(pn::NSNS) = pn.λ₂
