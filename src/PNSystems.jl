abstract type PNSystem{PNOrder,T} end
(PN::Type{PNS})(pnorder, t) where {PNS<:PNSystem} = PN{pnorder,t}()


mutable struct TaylorT1{PNOrder,T} <: PNSystem{PNOrder,T}
    M₁::T
    M₂::T
    χ⃗₁::QuatVec{T}
    χ⃗₂::QuatVec{T}
    R::Quaternion{T}
    v::T
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

"""
    recalculate!(u̇, u, pn)

Calculate the new values of `u̇` based on the values of `u`.  Note that this
modifies both `u̇` and `pn` in place.

"""
function recalculate!(u̇, u, pn::TaylorT1{PNOrder,T}) where {PNOrder,T}
    unpack!(pn, u)
    @unpack pn
    χ₁ = absvec(χ⃗₁)
    χ₂ = absvec(χ⃗₂)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(pn)
    let ℓ̂=ℓ̂(R), Ω⃗ᵪ₁=Ω⃗ᵪ₁(pn), Ω⃗ᵪ₂=Ω⃗ᵪ₂(pn), Ω⃗ₚ=Ω⃗ₚ(pn), 𝓕=𝓕(pn), 𝓔′=𝓔′(pn)
        v̇ = - (𝓕 + Ṁ₁ + Ṁ₂) / 𝓔′
        #v̇ = 2//5 * v^10 / (v/4)
        χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
        χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
        u̇[1] = Ṁ₁
        u̇[2] = Ṁ₂
        u̇[3:5] = ((Ṡ₁ / M₁^2 - 2χ₁ * Ṁ₁/M₁) * χ̂₁ + Ω⃗ᵪ₁ × χ⃗₁).vec
        u̇[6:8] = ((Ṡ₂ / M₂^2 - 2χ₂ * Ṁ₂/M₂) * χ̂₂ + Ω⃗ᵪ₂ × χ⃗₂).vec
        u̇[9:12] = (Ω⃗ₚ * R / 2).components
        u̇[13] = v̇
    end
    pn
end
