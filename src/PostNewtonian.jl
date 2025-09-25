module PostNewtonian

# Always explicitly address functions similar to functions defined in this package,
# which come from these packages:
using MacroTools: MacroTools
using FastDifferentiation: FastDifferentiation
using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions

# Otherwise, we just explicitly import specific functions:
using DataInterpolations: CubicSpline
using InteractiveUtils: methodswith
using LinearAlgebra: mul!
using Random: AbstractRNG, default_rng
using Quaternionic: QuatVec, Rotor, abs2vec, components, normalize, ⋅, ×
using SphericalFunctions: D!, Diterator, Dprep, Yiterator
using OrdinaryDiffEqVerner: Vern9
using SciMLBase:
    ODEFunction,
    ODEProblem,
    solve,
    remake,
    terminate!,
    CallbackSet,
    DiscreteCallback,
    VectorContinuousCallback,
    ODESolution,
    parameterless_type,
    FullSpecialize,
    AbstractDiffEqInterpolation,
    build_solution,
    get_du
using SciMLBase.ReturnCode: ReturnCode
using SymbolicIndexingInterface: SymbolCache
using RecursiveArrayTools: DiffEqArray
using StaticArrays: StaticArrays, SVector, MVector
using TestItems: @testitem

# See the "Code structure" section of the documentation for a description of the simple
# hierarchy into which this code is organized.  The different levels of that hierarchy are
# reflected cleanly in the files `include`d below.

# It's more common in PN to use `ln` — which I also prefer, as `log` seems ambiguous.
const ln = log

include("utilities.jl")
export termination_forwards,
    termination_backwards, dtmin_terminator, decreasing_v_terminator, nonfinite_terminator
using .MathConstants

include("pn_systems.jl")
export PNSystem, pn_order, BBH, BHBH, BHNS, NSNS, BNS, FDPNSystem, fd_pnsystem

include("pn_expansion.jl")
export PNExpansion, PNTerm, PNExpansionParameter

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
    n_hat,
    n̂,
    lambda_hat,
    λ̂,
    ell_hat,
    ℓ̂,
    Omega,
    Ω,
    S⃗₁,
    S⃗₂,
    S⃗,
    Σ⃗,
    χ⃗,
    χ⃗ₛ,
    χ⃗ₐ,
    chi_perp,
    χₚₑᵣₚ,
    chi_eff,
    χₑ,
    chi_p,
    χₚ,
    S⃗₀⁺,
    S⃗₀⁻,
    S₀⁺ₙ,
    S₀⁻ₙ,
    S₀⁺λ,
    S₀⁻λ,
    S₀⁺ₗ,
    S₀⁻ₗ,
    χ₁²,
    χ₂²,
    χ₁,
    χ₂,
    χ₁₂,
    χ₁ₗ,
    χ₂ₗ,
    χₛₗ,
    χₐₗ,
    Sₙ,
    Σₙ,
    Sλ,
    Σλ,
    Sₗ,
    Σₗ,
    sₗ,
    σₗ,
    S₁ₙ,
    S₁λ,
    S₁ₗ,
    S₂ₙ,
    S₂λ,
    S₂ₗ,
    rₕ₁,
    rₕ₂,
    Ωₕ₁,
    Ωₕ₂,
    sin²θ₁,
    sin²θ₂,
    ϕ̇̂₁,
    ϕ̇̂₂,
    Î₀₁,
    Î₀₂,
    κ₁,
    κ₂,
    κ₊,
    κ₋,
    λ₁,
    λ₂,
    λ₊,
    λ₋,
    Λ̃,
    Lambda_tilde

include("pn_expressions.jl")
export gw_energy_flux,
    𝓕,
    tidal_heating,
    binding_energy,
    𝓔,
    binding_energy_deriv,
    𝓔′,
    Omega_p,
    Ω⃗ₚ,
    Omega_chi1,
    Ω⃗ᵪ₁,
    Omega_chi2,
    Ω⃗ᵪ₂,
    #𝛡, aₗ, Ω⃗ᵪ  # Too obscure to bother with
    γₚₙ,
    inverse_separation,
    γₚₙ′,
    inverse_separation_deriv,
    γₚₙ⁻¹,
    inverse_separation_inverse,
    separation,  # r,
    separation_deriv,  # r′,
    separation_dot,  # ṙ,
    separation_inverse,  # r⁻¹,
    mode_weights!,
    h!,
    mode_weights_Ψ_M!,
    Ψ_M!

include("dynamics.jl")
export up_down_instability,
    estimated_time_to_merger, fISCO, ΩISCO, uniform_in_phase, orbital_evolution

include("waveforms.jl")
export coorbital_waveform,
    inertial_waveform,
    coorbital_waveform_computation_storage,
    inertial_waveform_computation_storage,
    coorbital_waveform!,
    inertial_waveform!

include("compatibility_layers.jl")
export GWFrames

include("assorted_binaries/examples.jl")
export superkick, hangup_kick
include("assorted_binaries/random.jl")
# Base.rand is the only function in that file, hence no need for exports

include("precompilation.jl")

include("predefinitions_Symbolics.jl")

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    # COV_EXCL_START

    function __init__()
        @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include(
            "../ext/PostNewtonianSymbolicsExt.jl"
        )
        @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" include(
            "../ext/PostNewtonianForwardDiffExt.jl"
        )
    end

    # COV_EXCL_STOP
end

end  # module PostNewtonian
