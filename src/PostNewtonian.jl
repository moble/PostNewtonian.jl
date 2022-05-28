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
export ℓ̂, n̂, λ̂, Ω_v, v_Ω

include("PNSystems.jl")
export PNSystem, TaylorT1

include("noneccentric_orbit.jl")
export noneccentric_evolution


"""
    PNParameters(M₁, M₂, S⃗₁, S⃗₂, R=Rotor(true), λ₁=false, λ₂=false)

Construct a `PNParameters` object, containing information about a PN system,
including masses and spins, orientation, and tidal-coupling constants.

"""
struct PNParameters{T<:Real}
    M₁::T
    M₂::T
    S⃗₁::QuatVec{T}
    S⃗₂::QuatVec{T}
    R::Rotor{T}
    λ₁::T
    λ₂::T
end
function PNParameters(M₁, M₂, S⃗₁::Real, S⃗₂::Real, R=Rotor(true), λ₁=false, λ₂=false)
    PNParameters(M₁, M₂, QuatVec(0,0,S⃗₁), QuatVec(0,0,S⃗₂), R, λ₁, λ₂)
end
function PNParameters(M₁, M₂, S⃗₁, S⃗₂, R=Rotor(true), λ₁=false, λ₂=false)
    T = promote_type(typeof(M₁), typeof(M₂), eltype(S⃗₁), eltype(S⃗₂), eltype(R), typeof(λ₁), typeof(λ₂))
    PNParameters{T}(T(M₁), T(M₂), QuatVec{T}(S⃗₁), QuatVec{T}(S⃗₂), Rotor{T}(R), T(λ₁), T(λ₂))
end



end
