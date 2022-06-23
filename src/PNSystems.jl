abstract type PNSystem{PNOrder,T} end
(PN::Type{PNS})(pnorder, t) where {PNS<:PNSystem} = PN{pnorder,t}()


mutable struct TaylorT1{PNOrder,T} <: PNSystem{PNOrder,T}
    TaylorT1{PNOrder,T}() where {PNOrder,T} = new()
end


"""
    recalculate!(u̇, u, pn)

Calculate the new values of `u̇` based on the values of `u`.  Note that this
modifies both `u̇` and `pn` in place.

"""
function recalculate!(u̇, u, pn::TaylorT1{PNOrder,T}) where {PNOrder,T}
    M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v = u
    χ⃗₁ = QuatVec(χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ)
    χ⃗₂ = QuatVec(χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ)
    R = Quaternion(Rʷ, Rˣ, Rʸ, Rᶻ)
    χ₁ = absvec(χ⃗₁)
    χ₂ = absvec(χ⃗₂)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(u)
    let ℓ̂=ℓ̂(R), Ω⃗ᵪ₁=Ω⃗ᵪ₁(u), Ω⃗ᵪ₂=Ω⃗ᵪ₂(u), Ω⃗ₚ=Ω⃗ₚ(u), 𝓕=𝓕(u), 𝓔′=𝓔′(u)
        Ω⃗ = Ω⃗ₚ + Ω(v=v, M=M₁+M₂) * ℓ̂
        v̇ = - (𝓕 + Ṁ₁ + Ṁ₂) / 𝓔′
        χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
        χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
        u̇[1] = Ṁ₁
        u̇[2] = Ṁ₂
        u̇[3:5] = ((Ṡ₁ / M₁^2 - 2χ₁ * Ṁ₁/M₁) * χ̂₁ + Ω⃗ᵪ₁ × χ⃗₁).vec
        u̇[6:8] = ((Ṡ₂ / M₂^2 - 2χ₂ * Ṁ₂/M₂) * χ̂₂ + Ω⃗ᵪ₂ × χ⃗₂).vec
        u̇[9:12] = (Ω⃗ * R / 2).components
        u̇[13] = v̇
    end
    pn
end
