module GWFrames

using ..PostNewtonian
using Quaternionic
using DataInterpolations: CubicSpline


"""
    PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; kwargs...)

Compute a PN waveform, with the same call signature as `GWFrames.PNWaveform`

This is essentially a compatibility layer for the corresponding function in the original
[`GWFrames`](https://github.com/moble/GWFrames/blob/01b39bfe/Code/PNWaveforms.cpp#L83-L88)
Python package, with several additional optional arguments: `inertial`, `dt`, `quiet`,
`ell_min`, `ell_max`, `Lambda1`, and `Lambda2` (see below).  Also, this function accepts
optional arguments either as positional arguments (which the original `GWFrames` requires)
or as keyword arguments.

!!! warning

    We do *not* expect the result of this function to be identical to the result from the
    `GWFrames` Python package.  In particular, this package uses more general expressions
    for the tidal-heating terms, fixes an error in the 2PN quadratic-spin terms for the
    waveform modes, uses more accurate (and efficient) ODE integration, and uses a more
    accurate method to compute the number of steps per orbit (by default).

    Also note that there are differences in the order of operations for computing
    intermediate variables.  Cancellation and roundoff error caused by these differences can
    have surprisingly large effects on the orbital evolution in many cases — particularly
    for precessing systems.  Results from the two packages have been painstakingly analyzed,
    leading to the conclusion that all differences are caused by such errors or the
    differences in formulations mentioned above.

The Julia interface is more detailed, flexible, and efficient than the simple `GWFrames`
interface that this function emulates.  In particular, [`orbital_evolution`](@ref) takes
essentially all the same arguments that
[`DifferentialEquations.solve`](https://diffeq.sciml.ai/dev/basics/common_solver_opts/)
takes, and returns a [`solution`](https://diffeq.sciml.ai/dev/basics/solution/) that
provides dense output and more details about the ODE solution itself.  For example, one
reason this function is more efficient than `GWFrames` is that we can use dense output to
solve with fewer timesteps, while accurately and efficiently interpolating to the requested
timesteps.  While `orbital_evolution` solves for the dynamics, [`coorbital_waveform`](@ref)
or [`inertial_waveform`](@ref) provides the actual waveform; both are returned by this
function.


## Required arguments

  * `Approximant`: Currently, only `"TaylorT1"` is supported.
  * `delta`: Fractional mass difference ``(M₁-M₂)/(M₁+M₂)``
  * `chi1_i`: Normalized spin vector ``S⃗₁/M₁²``
  * `chi2_i`: Normalized spin vector ``S⃗₂/M₂²``
  * `Omega_orb_i`: Orbital angular frequency at initial instant


## Optional arguments

As mentioned above, the following may be given *either* as positional arguments in this
order (though any number of them may be omitted from the end), or as keyword arguments.

  * `Omega_orb_0=Omega_orb_i`: Orbital angular frequency at first instant found in data.  If
    this is less than `Omega_orb_i`, the system is integrated backwards in time from the
    latter value to this value.
  * `R_frame_i=Rotor(1)`: Initial orientation of the frame.
  * `MinStepsPerOrbit=32`: Number of time steps in the output data per orbit.  Because the
    waveform modes go as high as ``m=8``, this number should be at least 16 to avoid Nyquist
    aliasing in those modes.  Note that this value may be overridden by `dt` (see below).
  * `PNWaveformModeOrder=4.0`: Maximum PN order of terms in the waveform formulas.
  * `PNOrbitalEvolutionOrder=4.0`: Maximum PN order of terms in the orbital-evolution
    formulas.
  * `inertial=false`: If `true`, transform waveform to the inertial frame; otherwise, the
    waveform will be in the co-orbital frame.
  * `dt=0`: Uniform time step size of the output.  If this is not a strictly positive
    number, `MinStepsPerOrbit` will be used instead.
  * `quiet=true`: If `false`, show informational messages about the reasons for terminating
    the ODE integration.  In either case, warnings will still be issued if terminating for
    bad or suspicious reasons.  See the documentation of [`orbital_evolution`](@ref) for an
    example of how to filter warnings also.
  * `ell_min=2`: The lowest ℓ value in the output waveform.
  * `ell_max=8`: The highest ℓ value in the output waveform.
  * `Lambda1=0`: Tidal-coupling parameter of object 1.
  * `Lambda2=0`: Tidal-coupling parameter of object 2.


## Returned values

This function returns a NamedTuple with the following keys:

  * `t`: The vector of time steps at which the data are evaluated.  The time ``t=0``
    corresponds to the initial values that are arguments to this function.
  * `data`: Matrix of complex values of the mode weights.  The shape is length(t) x 77.  The
    first dimension enumerates the values at each instant of time.  The second dimension
    enumerates the modes, starting with ``(2,-2)``, then ``(2,-1)``, up to ``(2,2)``,
    followed by ``(3,-3)``, and so on up to ``(8,8)``.  This is the opposite ordering as
    results from `GWFrames`, but the same as the ordering used by the `sxs` and `scri`
    packages.  However, also note that certain conversions between Julia and Python *may*
    transpose matrices, because Julia is Fortran-ordered by default, whereas numpy is
    C-ordered.  It is best to check the shape manually to be sure which dimension is which.
  * `frame`: Matrix of shape length(t) x 4 representing the frame-orientation quaternion as
    a function of time `t`.
  * `M1`, `M2`: Vectors of the respective masses as functions of time `t`.  Note that only
    at the time corresponding to `Omega_orb_i` will the total mass be precisely 1.
    Generally, tidal heating will lead to time-dependent masses.
  * `chi1`, `chi2`: Matrices of shape length(t) x 3 representing the spins as functions of
    time `t`.
  * `v`: PN velocity parameter as a function of time `t`.
  * `Phi`: Orbital phase as a function of time `t`.

Because this is a NamedTuple, the fields can be accessed much like the fields of a
`WaveformModes` object in the `scri` or `sxs` Python packages — as in `w.t` and `w.data`,
where `w` is the object returned by this function.
"""
function PNWaveform(
    Approximant::String, delta, chi1_i, chi2_i, Omega_orb_i,
    Omega_orb_0, R_frame_i=[1.0], MinStepsPerOrbit=32,
    PNWaveformModeOrder=4.0, PNOrbitalEvolutionOrder=4.0,
    inertial=false, dt=0.0, quiet=true,
    ell_min=2, ell_max=8, Lambda1=0, Lambda2=0,
    kwargs...
)
    # Note that this method's signature is missing the `Omega_orb_0` default
    # value; if it is not given, Julia selects the other (keyword-based)
    # method, which does have a default method for it.

    PNWaveform(
        Approximant, delta, chi1_i, chi2_i, Omega_orb_i;
        Omega_orb_0, R_frame_i, MinStepsPerOrbit,
        PNWaveformModeOrder, PNOrbitalEvolutionOrder,
        inertial, dt, quiet,
        ell_min, ell_max, Lambda1, Lambda2,
        kwargs...
    )
end

function PNWaveform(
    Approximant, delta, chi1_i, chi2_i, Omega_orb_i;
    Omega_orb_0=Omega_orb_i, R_frame_i=[1.0], MinStepsPerOrbit=32,
    PNWaveformModeOrder=4.0, PNOrbitalEvolutionOrder=4.0,
    inertial=false, dt=0.0, quiet=true,
    ell_min=2, ell_max=8, Lambda1=0, Lambda2=0,
    kwargs...
)
    M₁ = (1 + delta) / 2
    M₂ = (1 - delta) / 2
    χ⃗₁ = QuatVec(chi1_i)
    χ⃗₂ = QuatVec(chi2_i)
    Ωᵢ = Omega_orb_i
    Ω₁ = Omega_orb_0
    Rᵢ = Rotor(R_frame_i)

    saveat = if dt > 0.0
        if haskey(kwargs, :saveat)
            error("GWFrames.PNWaveform received both `dt` and `saveat` arguments")
        end
        dt
    elseif haskey(kwargs, :saveat)
        kwargs = Dict(kwargs)
        pop!(kwargs, :saveat)
    else
        []
    end

    # Inspiral
    solution = orbital_evolution(
        M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Λ₁=Lambda1, Λ₂=Lambda2,
        Ω₁=Ω₁, Rᵢ=Rᵢ,
        approximant=Approximant, PNOrder=PNOrbitalEvolutionOrder,
        quiet=quiet, saveat, kwargs...
    )
    if saveat == []
        solution = uniform_in_phase(solution, MinStepsPerOrbit)
    end

    # Waveform
    h = if inertial
        inertial_waveform(
            solution; ℓₘᵢₙ=ell_min, ℓₘₐₓ=ell_max, PNOrder=PNWaveformModeOrder
        )
    else
        coorbital_waveform(
            solution; ℓₘᵢₙ=ell_min, ℓₘₐₓ=ell_max, PNOrder=PNWaveformModeOrder
        )
    end

    # Return
    (
        t=solution.t,
        data=h',
        frame=solution[[:Rʷ, :Rˣ, :Rʸ, :Rᶻ], :]',
        M1=solution[:M₁],
        M2=solution[:M₂],
        chi1=solution[[:χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ], :]',
        chi2=solution[[:χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ], :]',
        v=solution[:v],
        Phi=solution[:Φ],
    )
end


end  # module GWFrames
