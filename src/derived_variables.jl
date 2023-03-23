module DerivedVariables

using ..PostNewtonian: PNSystem
using ..PostNewtonian.FundamentalVariables
using Quaternionic: 𝐢, 𝐣, 𝐤, QuatVec, (⋅), abs2vec

include("derived_variables/mass_combinations.jl")
export total_mass, M,
    reduced_mass, μ,
    reduced_mass_ratio, ν,
    mass_difference_ratio, δ,
    mass_ratio, q,
    chirp_mass, ℳ

include("derived_variables/orbital_elements.jl")
export n_hat, n̂,
    lambda_hat, λ̂,
    ell_hat, ℓ̂,
    Omega, Ω,
    lnv

include("derived_variables/spin_combinations.jl")
export S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ,
    χ₁², χ₂², χ₁, χ₂, χ₁₂,
    χ₁ₗ, χ₂ₗ, χₛₗ, χₐₗ,
    S₁ₙ, S₁λ, S₁ₗ, S₂ₙ, S₂λ, S₂ₗ

end
