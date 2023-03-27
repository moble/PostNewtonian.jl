module GWFrames

using ..PostNewtonian
using Quaternionic
using DataInterpolations: CubicSpline


"""
    PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; kwargs...)

Compute a PN waveform, with the same call signature as `GWFrames.PNWaveform`

This is essentially a compatibility layer for the corresponding function in the
original
[`GWFrames`](https://github.com/moble/GWFrames/blob/01b39bfe/Code/PNWaveforms.cpp#L83-L88)
Python package, with two additional optional arguments: `dt` and `quiet` (see
below).  Also, this function accepts optional arguments either as positional
arguments (which the original `GWFrames` requires) or as keyword arguments.

!!! warning

    We do *not* expect the result of this function to be identical to the
    result from the `GWFrames` Python package.  In particular, this package
    uses more general expressions for the tidal-heating terms, fixes an error
    in the 2PN quadratic-spin terms for the waveform modes, uses a more
    accurate method to compute the number of steps per orbit (by default), and
    uses more accurate (and efficient) ODE integration.

The Julia interface is more detailed, flexible, and efficient than the simple
`GWFrames` interface that this function emulates.  In particular,
[`inspiral`](@ref) takes essentially all the same arguments that
[`DifferentialEquations.solve`](https://diffeq.sciml.ai/dev/basics/common_solver_opts/)
takes, and returns a [`solution`](https://diffeq.sciml.ai/dev/basics/solution/)
that provides dense output and more details about the ODE solution itself.  For
example, one reason this function is more efficient than `GWFrames` is that we
can use dense output to solve with fewer timesteps, while accurately and
efficiently interpolating to the requested timesteps.  While `inspiral` solves
for the dynamics, [`h!`](@ref) provides the actual mode weights, which are also
returned by this function.


## Required arguments

  * `Approximant`: Currently, only `"TaylorT1"` is supported.
  * `delta`: Fractional mass difference ``(M₁-M₂)/(M₁+M₂)``
  * `chi1_i`: Normalized spin vector ``S⃗₁/M₁²``
  * `chi2_i`: Normalized spin vector ``S⃗₂/M₂²``
  * `Omega_orb_i`: Orbital angular frequency at initial instant


## Optional arguments

As mentioned above, the following may be given *either* as positional arguments
in this order (though any number of them may be omitted from the end), or as
keyword arguments.

  * `Omega_orb_0=Omega_orb_i`: Orbital angular frequency at first instant found
    in data.  If this is less than `Omega_orb_i`, the system is integrated
    backwards in time from the latter value to this value.
  * `R_frame_i=Rotor(1)`: Initial orientation of the frame.
  * `MinStepsPerOrbit=32`: Number of time steps in the output data per orbit.
    Because the waveform modes go as high as ``m=8``, this number should be at
    least 16 to avoid Nyquist aliasing in those modes.  Note that this value
    may be overridden by `dt` (see below).
  * `PNWaveformModeOrder=3.5`: Maximum PN order of terms in the
    waveform formulas.  Currently, only 3.5 is supported.
  * `PNOrbitalEvolutionOrder=4.0`: Maximum PN order of terms in the
    orbital-evolution formulas.  Currently, only 4.0 is supported.
  * `dt=0`: Uniform time step size of the output.  If this is not a strictly
    positive number, `MinStepsPerOrbit` will be used instead.
  * `quiet=true`: If `false`, show informational messages about the reasons for
    terminating the ODE integration.  In either case, warnings will still be
    issued if terminating for bad or suspicious reasons.  See the documentation
    of [`inspiral`](@ref) for an example of how to filter warnings also.


## Returned values

This function returns a NamedTuple with the following keys:

  * `t`: The vector of time steps at which the data are evaluated.  The time
    ``t=0`` corresponds to the initial values that are arguments to this
    function.
  * `data`: Matrix of complex values of the mode weights.  The shape is
    length(t) x 77.  The first dimension enumerates the values at each instant
    of time.  The second dimension enumerates the modes, starting with
    ``(2,-2)``, then ``(2,-1)``, up to ``(2,2)``, followed by ``(3,-3)``, and
    so on up to ``(8,8)``.  This is the opposite ordering as results from
    `GWFrames`, but the same as the ordering used by the `sxs` and `scri`
    packages.  However, also note that certain conversions between Julia and
    Python *may* transpose matrices, because Julia is Fortran-ordered by
    default, whereas numpy is C-ordered.  It is best to check the shape
    manually to be sure which dimension is which.
  * `frame`: Matrix of shape length(t) x 4 representing the frame-orientation
    quaternion as a function of time `t`.
  * `M1`, `M2`: Vectors of the respective masses as functions of time `t`.
    Note that only at the time corresponding to `Omega_orb_i` will the total
    mass be precisely 1.  Generally, tidal heating will lead to time-dependent
    masses.
  * `chi1`, `chi2`: Matrices of shape length(t) x 3 representing the spins as
    functions of time `t`.
  * `v`: PN velocity parameter as a function of time `t`.
  * `Phi`: Orbital phase as a function of time `t`.

Because this is a NamedTuple, the fields can be accessed much like the fields
of a `WaveformModes` object in the `scri` or `sxs` Python packages — as in
`w.t` and `w.data`, where `w` is the object returned by this function.

"""
function PNWaveform(
    Approximant::String, delta, chi1_i, chi2_i, Omega_orb_i,
    Omega_orb_0, R_frame_i=[1.0], MinStepsPerOrbit=32,
    PNWaveformModeOrder=3.5, PNOrbitalEvolutionOrder=4.0, dt=0.0, quiet=true
)
    # Note that this method's signature is missing the `Omega_orb_0` default
    # value; if it is not given, Julia selects the other (keyword-based)
    # method, which does have a default method for it.

    PNWaveform(
        Approximant, delta, chi1_i, chi2_i, Omega_orb_i;
        Omega_orb_0, R_frame_i, MinStepsPerOrbit,
        PNWaveformModeOrder, PNOrbitalEvolutionOrder, dt, quiet
    )
end
function PNWaveform(
    Approximant, delta, chi1_i, chi2_i, Omega_orb_i;
    Omega_orb_0=Omega_orb_i, R_frame_i=[1.0], MinStepsPerOrbit=32,
    PNWaveformModeOrder=3.5, PNOrbitalEvolutionOrder=4.0, dt=0.0, quiet=true,
    ell_min=2, ell_max=8, lambda1=0, lambda2=0
)
    if PNWaveformModeOrder ≉ 3.5
        @error "`PNWaveformModeOrder` other than 3.5 is not yet supported"
    end
    M₁ = (1 + delta) / 2
    M₂ = (1 - delta) / 2
    χ⃗₁ = QuatVec(chi1_i)
    χ⃗₂ = QuatVec(chi2_i)
    Ωᵢ = Omega_orb_i
    Ω₁ = Omega_orb_0
    Rᵢ = Rotor(R_frame_i)

    # Inspiral
    solution = inspiral(
        M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; λ₁=lambda1, λ₂=lambda2,
        Ω₁=Ω₁, Rᵢ=Rᵢ,
        approximant=Approximant, PNOrder=PNOrbitalEvolutionOrder,
        quiet=quiet,
        saveat=dt > 0 ? dt : []
    )
    if dt ≤ 0
        solution = let
            Φ = solution[14, :]
            t = solution.t
            δΦ = 2π / MinStepsPerOrbit
            Φrange = range(extrema(Φ)..., step=δΦ)
            t_Φ = CubicSpline(t, Φ)(Φrange)
            solution(t_Φ)
        end
    end

    # Modes
    ℓmin = ell_min
    ℓmax = ell_max
    n_modes = (ℓmax+1)^2 - ℓmin^2
    pnsystem = if iszero(lambda1) && iszero(lambda2)
        BBH(copy(solution.u[1]))
    elseif iszero(lambda1)
        BHNS(copy(solution.u[1]), lambda2)
    else
        NSNS(copy(solution.u[1]), lambda1, lambda2)
    end
    h = Matrix{Complex{eltype(solution)}}(undef, n_modes, length(solution.t))
    for (hi, ui) in zip(axes(h, 2), solution.u)
        pnsystem.state[:] .= ui
        h!(@view(h[:, hi]), pnsystem; ℓmin, ℓmax)
    end

    # Return
    (
        t=solution.t,
        data=h',
        frame=solution[9:12, :]', # R
        M1=solution[1, :], # M₁
        M2=solution[2, :], # M₂
        chi1=solution[3:5, :]', # χ⃗₁
        chi2=solution[6:8, :]', # χ⃗₂
        v=solution[13, :], # v
        Phi=solution[14, :], # Φ
    )
end


end  # module GWFrames
