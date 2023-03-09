# abstract type PNSystem{PNOrder,T} end
# (PN::Type{PNS})(pnorder, t) where {PNS<:PNSystem} = PN{pnorder,t}()


# mutable struct TaylorT1{PNOrder,T} <: PNSystem{PNOrder,T}
#     TaylorT1{PNOrder,T}() where {PNOrder,T} = new()
# end


"""
    recalculate!(u̇, u, p)

Calculate the new values of `u̇` based on the values of `u`.

"""
function recalculate!(u̇, u, p)
    M₁, M₂, χ⃗₁, χ⃗₂, R, v = (
        u[1], u[2], QuatVec(u[3:5]...), QuatVec(u[6:8]...), Rotor(u[9:12]...), u[13]
    )
    χ₁, χ₂ = absvec(χ⃗₁), absvec(χ⃗₂)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(u)
    let ℓ̂=ℓ̂(R), Ω=Ω(v=v, M=M₁+M₂), Ω⃗ₚ=Ω⃗ₚ(u),
        Ω⃗ᵪ₁=Ω⃗ᵪ₁(u), Ω⃗ᵪ₂=Ω⃗ᵪ₂(u), 𝓕=𝓕(u), 𝓔′=𝓔′(u)
        Ω⃗ = Ω⃗ₚ + Ω * ℓ̂
        v̇ = - (𝓕 + Ṁ₁ + Ṁ₂) / 𝓔′
        χ̂₁ = ifelse(iszero(χ₁), ℓ̂, χ⃗₁ / χ₁)
        χ̂₂ = ifelse(iszero(χ₂), ℓ̂, χ⃗₂ / χ₂)
        u̇[1] = Ṁ₁
        u̇[2] = Ṁ₂
        u̇[3:5] = vec((Ṡ₁ / M₁^2 - 2χ₁ * Ṁ₁/M₁) * χ̂₁ + Ω⃗ᵪ₁ × χ⃗₁)
        u̇[6:8] = vec((Ṡ₂ / M₂^2 - 2χ₂ * Ṁ₂/M₂) * χ̂₂ + Ω⃗ᵪ₂ × χ⃗₂)
        u̇[9:12] = components(Ω⃗ * R / 2)
        u̇[13] = v̇
        if length(u̇) == 14
            u̇[14] = Ω
        end
    end
    nothing
end
