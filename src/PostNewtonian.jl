module PostNewtonian

# Always explicitly address functions similar to functions defined in this package,
# which come from these packages:
import MacroTools
import Symbolics
import SymbolicUtils
import FastDifferentiation

# Otherwise, we just explicitly import specific functions:
using DataInterpolations: CubicSpline
using InteractiveUtils: methodswith
using LinearAlgebra: mul!
using OrdinaryDiffEq: Vern9, #AutoVern9, Rodas5P,
    ODEFunction, ODEProblem, solve, remake,
    terminate!, CallbackSet, DiscreteCallback, VectorContinuousCallback
using Quaternionic: QuatVec, Rotor, abs2vec, components, normalize, (⋅), (×)
using Random: AbstractRNG, GLOBAL_RNG
using RecursiveArrayTools: DiffEqArray
using SciMLBase: ODESolution, parameterless_type, FullSpecialize,
    AbstractDiffEqInterpolation, build_solution, get_du
using SciMLBase.ReturnCode: ReturnCode
using SphericalFunctions: D!, Diterator, Dprep, Yiterator
using SymbolicIndexingInterface: SymbolCache
using RuntimeGeneratedFunctions: get_expression

# See the "Code structure" section of the documentation for a description of the simple
# hierarchy into which this code is organized.  The different levels of that hierarchy are
# reflected cleanly in the files `include`d below.


include("utilities.jl")
export termination_forwards, termination_backwards,
    dtmin_terminator, decreasing_v_terminator, nonfinite_terminator
using .MathConstants


include("systems.jl")
export PNSystem, BBH, BHNS, NSNS, SymbolicPNSystem, symbolic_pnsystem, FDPNSystem, fd_pnsystem, pn_order


include("fundamental_variables.jl")
using .FundamentalVariables
#export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₁, Λ₂  # Avoid clashes: don't export


include("derived_variables.jl")
using .DerivedVariables
export total_mass,  # M,  # Avoid clashes: don't export nicer names for important variables
    reduced_mass,  # μ,
    reduced_mass_ratio,  # ν,
    mass_difference_ratio,  # δ,
    mass_ratio,  # q,
    chirp_mass,  # ℳ,
    # X1, X₁,
    # X2, X₂,
    n_hat, n̂,
    lambda_hat, λ̂,
    ell_hat, ℓ̂,
    Omega, Ω,
    S⃗₁, S⃗₂, S⃗, Σ⃗, χ⃗, χ⃗ₛ, χ⃗ₐ,
    chi_perp, χₚₑᵣₚ,
    chi_eff, χₑ,
    chi_p, χₚ,
    Sₙ, Σₙ, Sλ, Σλ, Sₗ, Σₗ


include("pn_expressions.jl")
export gw_energy_flux, 𝓕,
    tidal_heating,
    binding_energy, 𝓔,
    binding_energy_deriv, 𝓔′,
    Omega_p, Ω⃗ₚ,
    Omega_chi1, Ω⃗ᵪ₁,
    Omega_chi2, Ω⃗ᵪ₂,
    #𝛡, γₚₙ, aₗ, Ω⃗ᵪ  # Too obscure to bother with
    mode_weights!, h!


include("dynamics.jl")
export up_down_instability, estimated_time_to_merger, fISCO, ΩISCO,
    uniform_in_phase, orbital_evolution


include("waveforms.jl")
export coorbital_waveform, inertial_waveform,
    coorbital_waveform_computation_storage, inertial_waveform_computation_storage,
    coorbital_waveform!, inertial_waveform!


include("compatibility_layers.jl")
export GWFrames


include("assorted_binaries/examples.jl")
export superkick, hangup_kick
include("assorted_binaries/random.jl")
# Base.rand is the only function in that file


include("precompilation.jl")


end  # module PostNewtonian
