module PostNewtonian

using StaticArrays
using Quaternionic
using Symbolics
using SymbolicUtils
using OrdinaryDiffEq

include("constants.jl")

include("PNSystems.jl")
export PNSystem, TaylorT1

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

include("pn_dynamics/tidal_heating.jl")
export tidal_heating

include("pn_dynamics/precession.jl")
export Ω⃗ₚ, Omega_p,
    Ω⃗ᵪ₁, Omega_chi1,
    Ω⃗ᵪ₂, Omega_chi2,
    Ω⃗ᵪ, 𝛡, γ, aₗ

include("pn_dynamics/flux.jl")
export 𝓕, gravitational_wave_flux,
    𝓕EMRI, gravitational_wave_flux_EMRI,
    𝓕NS, gravitational_wave_flux_NS

include("pn_dynamics/binding_energy.jl")
export 𝓔, binding_energy,
    𝓔′, binding_energy_deriv,
    𝓔NS, binding_energy_NS

include("up_down_instability.jl")
export up_down_instability

include("noneccentric_orbit.jl")
export noneccentric_evolution




end
