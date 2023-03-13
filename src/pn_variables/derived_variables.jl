module DerivedVariables

using ..PostNewtonian: PNSystem
using Quaternionic: 𝐢, 𝐣, 𝐤

include("mass_combinations.jl")
export M, total_mass,
    μ, reduced_mass,
    ν, reduced_mass_ratio,
    δ, mass_difference_ratio,
    q, mass_ratio,
    ℳ, chirp_mass

include("orbital_elements.jl")
export n̂, λ̂, ℓ̂, Ω

include("spin_combinations.jl")
export S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ

end
