"""
    recalculate!(u̇, u, p)

Calculate the new values of `u̇` based on the values of `u`.

"""
@compute_pn_variables 3 function recalculate!(u̇, u, p)
    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂) = tidal_heating(p)
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
    if length(u̇) > 13
        u̇[14] = Ω
    end
    nothing
end
