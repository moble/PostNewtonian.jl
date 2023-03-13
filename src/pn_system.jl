abstract type PNSystem{T, PNOrder} end

eltype(::PNSystem{T, PNOrder}) where {T, PNOrder} = T
pn_order(::PNSystem{T, PNOrder}) where {T, PNOrder} = PNOrder

struct TaylorT1{T, PNOrder} <: PNSystem{T, PNOrder}
    u::AbstractVector{T}
end
function TaylorT1(;M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing, PNOrder=typemax(Int))
    u = [M₁; M₂; vec(χ⃗₁); vec(χ⃗₂); components(R); v]
    if !isnothing(Φ)
        u = [u; Φ]
    end
    T = eltype(u)
    TaylorT1{T, Rational(PNOrder)}(u)
end
