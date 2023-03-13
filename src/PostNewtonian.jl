module PostNewtonian

using StaticArrays
using Quaternionic
using MacroTools
using Symbolics
using SymbolicUtils
using SciMLBase
using DiffEqBase
using RecursiveArrayTools
using OrdinaryDiffEq

include("utilities/misc.jl")

include("utilities/mathconstants.jl")
using .MathConstants

include("pn_systems/pn_system.jl")
export AbstractPNSystem, BBH

include("pn_variables/fundamental_variables.jl")
using .FundamentalVariables
#export M₁, M₂, χ⃗₁, χ⃗₂, R, v

include("pn_variables/derived_variables.jl")
using .DerivedVariables
export total_mass, reduced_mass, reduced_mass_ratio,
    mass_difference_ratio, mass_ratio, chirp_mass,
    n̂, λ̂, ℓ̂, Ω,
    S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ

include("utilities/macros.jl")

include("pn_expressions/tidal_heating.jl")
export tidal_heating

include("pn_expressions/precession.jl")
export Ω⃗ₚ, Omega_p,
    Ω⃗ᵪ₁, Omega_chi1,
    Ω⃗ᵪ₂, Omega_chi2,
    𝛡, γ, aₗ

include("pn_expressions/flux.jl")
export 𝓕, gw_energy_flux

include("pn_expressions/binding_energy.jl")
export 𝓔, binding_energy,
    𝓔′, binding_energy_deriv

include("pn_expressions/mode_weights.jl")
export h!, mode_weights!

include("pn_dynamics/up_down_instability.jl")
export up_down_instability

include("utilities/termination_criteria.jl")
export termination_forwards, termination_backwards,
    dtmin_terminator, nonfinite_terminator

include("pn_dynamics/PNSystems.jl")

include("utilities/combine_solutions.jl")

include("pn_dynamics/inspiral.jl")
export inspiral

include("compatibility_layers/gwframes.jl")
export GWFrames

end  # module PostNewtonian
