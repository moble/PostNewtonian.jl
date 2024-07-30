"""
    superkick(;v=0.2, χ=0.99, PNOrder=typemax(Int))
    superkick(;v=0.2, chi=0.99, PNOrder=typemax(Int))

Construct a black-hole binary in "superkick" configuration.

This is the scenario first published by [Campanelli et al.
(2007)](https://arxiv.org/abs/gr-qc/0701164), which has equal-mass black holes with spins of
equal magnitude oriented in opposite directions in the orbital plane.  This configuration
produces large asymmetrical emission of gravitational-wave linear momentum along the ``+z``
or ``-z`` directions, depending on which part of the orbit the binary is in.  Depending on
when the system mergers, the remnant may then acquire a huge recoil velocity.

(That recoil velocity, of course, depends on details of the merger, which post-Newtonian
approximations cannot describe.  Therefore, we cannot actually use PN methods to derive any
useful information about the remnant.)

Note that the name "superkick" appeared in the literature before a class of systems that can
achieve even larger recoil velocities: see [`hangup_kick`](@ref) for such systems.

# Examples
```julia-repl
julia> pnsystem = superkick(v=0.1)
BBH{Vector{Float64}, 9223372036854775805//2}([0.5, 0.5, 0.99, 0.0, 0.0, -0.99, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.1, 0.0])

julia> inspiral = orbital_evolution(pnsystem);

julia> inspiral[:v, 1]
0.1

julia> absvec(PostNewtonian.χ⃗₁(inspiral[1]))
0.99
```
"""
function superkick(; v=0.2, chi=0.99, χ=chi, PNOrder=typemax(Int))
    M₁ = M₂ = 0.5
    χ⃗₁ = [χ, 0, 0]
    χ⃗₂ = [-χ, 0, 0]
    R = Rotor(1)
    return BBH(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, PNOrder)
end

@doc raw"""
    hangup_kick(;v=0.2, χ=0.99, θ=deg2rad(50.98), ϕ=deg2rad(30.0), PNOrder=typemax(Int))
    hangup_kick(;v=0.2, chi=0.99, theta=deg2rad(50.98), phi=deg2rad(30.0), PNOrder=typemax(Int))

Construct a black-hole binary in [hangup-kick](https://arxiv.org/abs/1908.04382)
configuration.

The spin magnitudes are both equal to `χ`.  The direction of ``\vec{\chi}_1`` is given by
the spherical coordinates (`θ`, `ϕ`).  ``\vec{\chi}_2`` is the same, except that its
projection into the orbital plane is opposite to that of ``\vec{\chi}_1``.

See also [`superkick`](@ref) for the original systems, which can't actually achieve recoil
kicks as large as those achieved by "hangup-kick" systems.

The physical mechanism here is similar to the one described in [`superkick`](@ref), with an
additional "hangup" due to the fact that the spins have components parallel to the orbital
angular velocity.  The physical interpretation of that part is that any system with such
spin components will not merge as quickly because the remnant black hole must have total
dimensionless spin less than 1, so excess angular momentum must be radiated away.  This
gives the hangup-kick configuration more time to develop a large recoil.

# Examples
```julia-repl
julia> pnsystem = hangup_kick(v=0.1)
BBH{Vector{Float64}, 9223372036854775805//2}([0.5, 0.5, 0.3845784887294712, 0.6661094819774992, 0.6232957115416596, -0.3845784887294712, -0.6661094819774992, 0.6232957115416596, 1.0, 0.0, 0.0, 0.0, 0.1, 0.0])

julia> inspiral = orbital_evolution(pnsystem);

julia> inspiral[:v, 1]
0.1

julia> absvec(PostNewtonian.χ⃗₁(inspiral[1]))
0.99
```
"""
function hangup_kick(;
    v=0.2,
    chi=0.99,
    χ=chi,
    theta=deg2rad(50.98),
    θ=theta,
    phi=deg2rad(30.0),
    ϕ=phi,
    PNOrder=typemax(Int),
)
    M₁ = M₂ = 0.5
    χ⃗₁ = χ * [sin(ϕ) * sin(θ), cos(ϕ) * sin(θ), cos(θ)]
    χ⃗₂ = χ * [-sin(ϕ) * sin(θ), -cos(ϕ) * sin(θ), cos(θ)]
    R = Rotor(1)
    return BBH(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, PNOrder)
end
