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
    χ⃗̇₁::QuatVec{T}
    χ⃗̇₂::QuatVec{T}
    Ṙ::Quaternion{T}
    v̇::T
    TaylorT1{PNOrder,T}() where {PNOrder,T} = new()
end

macro unpack(pn)
    esc(quote
            M₁ = $(pn).M₁
            M₂ = $(pn).M₂
            χ⃗₁ = $(pn).χ⃗₁
            χ⃗₂ = $(pn).χ⃗₂
            R = $(pn).R
            v = $(pn).v
    end)
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
    @unpack pn
    χ₁ = absvec(χ⃗₁)
    χ₂ = absvec(χ⃗₂)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(pn)
    let ℓ̂=ℓ̂(R), Ω=Ω(v=v, M=M₁+M₂)
        χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
        χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
        pn.Ṁ₁ = Ṁ₁
        pn.Ṁ₂ = Ṁ₂
        pn.χ⃗̇₁ = (Ṡ₁ / M₁^2) * χ̂₁
        pn.χ⃗̇₂ = (Ṡ₂ / M₂^2) * χ̂₂
        pn.Ṙ = Ω * ℓ̂ * R / 2
        pn.v̇ = √(v)
    end
    pn
end
