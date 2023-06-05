module DerivedVariables

using ..PostNewtonian: PNSystem, VecOrPNSystem
using ..PostNewtonian.FundamentalVariables
using Quaternionic: 𝐢, 𝐣, 𝐤, QuatVec, (⋅), (×), abs2vec, absvec

include("derived_variables/mass_combinations.jl")
export total_mass, M,
    reduced_mass, μ,
    reduced_mass_ratio, ν,
    mass_difference_ratio, δ,
    mass_ratio, q,
    chirp_mass, ℳ,
    X1, X₁,
    X2, X₂

include("derived_variables/orbital_elements.jl")
export n_hat, n̂,
    lambda_hat, λ̂,
    ell_hat, ℓ̂,
    Omega, Ω,
    lnv

include("derived_variables/spin_combinations.jl")
export S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    chi_perp, χₚₑᵣₚ, χₑ, chi_eff, χₚ, chi_p,
    S⃗₀⁺, S⃗₀⁻, S₀⁺ₗ, S₀⁻ₗ,
    χ₁², χ₂², χ₁, χ₂, χ₁₂,
    χ₁ₗ, χ₂ₗ, χₛₗ, χₐₗ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ, sₗ, σₗ,
    S₁ₙ, S₁λ, S₁ₗ, S₂ₙ, S₂λ, S₂ₗ

include("derived_variables/horizons.jl")
export rₕ₁, rₕ₂, Ωₕ₁, Ωₕ₂,
    sin²θ₁, sin²θ₂, ϕ̇̂₁, ϕ̇̂₂,
    Î₀₁, Î₀₂,
    κ₁, κ₂, κ₊, κ₋,
    λ₁, λ₂, λ₊, λ₋

end
