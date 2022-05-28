abstract type PNSystem end

mutable struct TaylorT1{T} <: PNSystem
    ṁ₁::T
    ṁ₂::T
    Ω⃗ᵪ₁::QuatVec{T}
    Ω⃗ᵪ₂::QuatVec{T}
    Ω⃗::QuatVec{T}
    v̇::T
end
TaylorT1{T}(u) = TaylorT1{T}(u[1], u[2], u[3:5], u[6:8], u[9:11], u[12])

function recalculate!(pn::TaylorT1{T}, u) where T
    pn.ṁ₁ = 0
    pn.ṁ₂ = 0
    pn.Ω⃗ᵪ₁ = QuatVec{T}(0)
    pn.Ω⃗ᵪ₂ = QuatVec{T}(0)
    pn.Ω⃗ = QuatVec(0, 0, Ω_v(u[end], u[1]+u[2]))
    pn.v̇ = u[end]
    nothing
end
