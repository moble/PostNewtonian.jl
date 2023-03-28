module PostNewtonian

using Random: AbstractRNG, GLOBAL_RNG
using InteractiveUtils: methodswith
using StaticArrays
using Quaternionic
using RecursiveArrayTools
using MacroTools
using Symbolics
using SymbolicUtils
using SciMLBase
using DiffEqBase
using OrdinaryDiffEq

# See the "Code structure" section of the documentation for a description of the simple
# hierarchy into which this code is organized.  The different levels of that hierarchy are
# reflected cleanly in the files `include`d below.


include("utilities.jl")
export termination_forwards, termination_backwards,
    dtmin_terminator, decreasing_v_terminator, nonfinite_terminator
using .MathConstants


include("systems.jl")
export PNSystem, BBH, BHNS, NSNS, SymbolicPNSystem, symbolic_pnsystem, pn_order


include("fundamental_variables.jl")
using .FundamentalVariables
#export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, λ₁, λ₂  # Avoid clashes: don't export


include("derived_variables.jl")
using .DerivedVariables
export total_mass,  # M,  # Avoid clashes: don't export nicer names for important variables
    reduced_mass,  # μ,
    reduced_mass_ratio,  # ν,
    mass_difference_ratio,  # δ,
    mass_ratio,  # q,
    chirp_mass,  # ℳ,
    n_hat, n̂,
    lambda_hat, λ̂,
    ell_hat, ℓ̂,
    Omega, Ω,
    S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ


include("pn_expressions.jl")
export gw_energy_flux, 𝓕,
    tidal_heating,
    binding_energy, 𝓔,
    binding_energy_deriv, 𝓔′,
    Omega_p, Ω⃗ₚ,
    Omega_chi1, Ω⃗ᵪ₁,
    Omega_chi2, Ω⃗ᵪ₂,
    #𝛡, γₚ, aₗ, Ω⃗ᵪ  # Too obscure to bother with
    mode_weights!, h!


include("dynamics.jl")
export estimated_time_to_merger, up_down_instability, orbital_evolution


include("evaluation.jl")


include("compatibility_layers.jl")
export GWFrames


include("examples.jl")


end  # module PostNewtonian
