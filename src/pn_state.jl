abstract type PNState{T, PNOrder} end

pnorder(::PNState{T, PNOrder}) where {T, PNOrder} = PNOrder
eltype(::PNState{T, PNOrder}) where {T, PNOrder} = T

struct TaylorT1{T, PNOrder} <: PNState{T, PNOrder}
    u::AbstractVector{T}
end
