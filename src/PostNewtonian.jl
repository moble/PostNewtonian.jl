module PostNewtonian

using StaticArrays
using Quaternionic
using Symbolics
using DifferentialEquations


include("masses.jl")
export μ, reduced_mass,
    ν, reduced_mass_ratio,
    δ, mass_difference_ratio,
    q, mass_ratio,
    ℳ, chirp_mass

include("spins.jl")
export χ⃗, S, Σ, χₛ, χₐ

include("orbital_elements.jl")
export ℓ̂, n̂, λ̂, Ω, v

include("up_down_instability.jl")
export up_down_instability

include("PNSystems.jl")
export PNSystem, TaylorT1

include("pn_dynamics/tidal_heating.jl")
export tidal_heating

include("pn_dynamics/precession.jl")
export Ω⃗ₚ, Ω⃗ᵪ₁, Ω⃗ᵪ₂, Ω⃗ᵪ, 𝛡, γ, aₗ

include("noneccentric_orbit.jl")
export noneccentric_evolution


end
