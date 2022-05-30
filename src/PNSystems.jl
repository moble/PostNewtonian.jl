abstract type PNSystem{PNOrder,T} end
(PN::Type{PNS})(pnorder, t) where {PNS<:PNSystem} = PN{pnorder,t}()


mutable struct TaylorT1{PNOrder,T} <: PNSystem{PNOrder,T}
    M₁::T
    M₂::T
    χ⃗₁::QuatVec{T}
    χ⃗₂::QuatVec{T}
    R::Quaternion{T}
    v::T
    Ṁ₁::T
    Ṁ₂::T
    Ω⃗ᵪ₁::QuatVec{T}
    Ω⃗ᵪ₂::QuatVec{T}
    Ω⃗::Quaternion{T}
    v̇::T
    TaylorT1{PNOrder,T}() where {PNOrder,T} = new()
end

function unpack!(pn::PNSystem{PNOrder,T}, u) where {PNOrder,T}
    pn.M₁ = u[1]
    pn.M₂ = u[2]
    pn.χ⃗₁ = QuatVec{T}(u[3:5]...)
    pn.χ⃗₂ = QuatVec{T}(u[6:8]...)
    pn.R = Quaternion{T}(u[9:12]...)
    pn.v = u[13]
    pn
end

function recalculate!(pn::TaylorT1{PNOrder,T}, u) where {PNOrder,T}
    unpack!(pn, u)
    pn.Ṁ₁ = 0
    pn.Ṁ₂ = 0
    pn.Ω⃗ᵪ₁ = QuatVec{T}(0)
    pn.Ω⃗ᵪ₂ = QuatVec{T}(0)
    pn.Ω⃗ = Ω(v=pn.v, M=pn.M₁+pn.M₂) * ℓ̂(pn.R)
    pn.v̇ = √(pn.v)
    pn
end
