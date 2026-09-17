"""
    estimated_time_to_merger(M, ν, v)
    estimated_time_to_merger(pnsystem)

Compute the lowest-order PN approximation for the time to merger, starting from PN velocity
parameter `v`.

This is used internally as a convenient way to estimate how long the inspiral integration
should run for; we don't want it to integrate forever if PN has broken down.  However, it
can be a very poor approximation, especially close to merger, and doubly so if the spins or
eccentricity are significant.
"""
function estimated_time_to_merger(M, ν, v)
    return 5M / (256ν * v^8)
end

function estimated_time_to_merger(pnsystem)
    return estimated_time_to_merger(M(pnsystem), ν(pnsystem), v(pnsystem))
end

"""
    fISCO(q, M)
    fISCO(pnsystem)

Compute the "BKL" approximation for the ISCO (Innermost Stable Circular Orbit) frequency.

This is taken from Eq. (5) of [Hanna et al. (2008)](https://arxiv.org/abs/0801.4297).  Note
that this does not account for the spins of the objects in the binary, so that this returns
a very crude estimate of a frequency of interest.
"""
function fISCO(q, M)
    let π = oftype(q, π)
        if q > 1
            q = 1 / q
        end
        (10 + q * (28 + q * (-26 + q * 8))) / (10π * (6M)^(3//2))
    end
end
function fISCO(pnsystem)
    return fISCO(q(pnsystem), M(pnsystem))
end

"""
    ΩISCO(q,M)
    ΩISCO(pnsystem)

2π times [`fISCO`](@ref).
"""
function ΩISCO(q, M)
    return 2oftype(q, π) * fISCO(q, M)
end
function ΩISCO(pnsystem)
    return 2eltype(pnsystem)(π) * fISCO(pnsystem)
end

@doc raw"""
    uniform_in_phase(solution, saves_per_orbit)

Interpolate `solution` to uniform steps in phase.

By default, the `solution` returned by [`orbital_evolution`](@ref) may be sampled very
sparsely — too sparsely to satisfy the Nyquist limit of the waveform.  If the waveform
extends to ``\ell_{\mathrm{max}}``, there will be modes varying slightly more rapidly than
``\exp\left(\pm i\, \ell_{\mathrm{max}}\, \Phi \right)``, where ``\Phi`` is the orbital
phase.  If the frequency were constant, this would require at least ``2\ell_{\mathrm{max}}``
samples per orbit.  To incorporate a safety factor, ``4\ell_{\mathrm{max}}`` seems to work
fairly reliably.

See also the `saves_per_orbit` and `saveat` arguments to [`orbital_evolution`](@ref), as
well as interpolation-in-time capabilities of the result of that function.
"""
function uniform_in_phase(solution, saves_per_orbit)
    let π = eltype(solution)(π)
        t = solution.t
        PNType = typeof(solution.prob.p)
        Φ = solution[:Φ, :]
        δΦ = 2π / saves_per_orbit
        Φrange = range(extrema(Φ)...; step=δΦ)
        t_Φ = CubicSpline(t, Φ)(Φrange)
        # Ensure that t=0 is interpolated back
        # to *exactly* t=0 instead of, e.g., -1e-24:
        t_Φ[1] = t[1]
        solution(t_Φ)
    end
end

function default_termination_criteria_forwards(pnsystem, vₑ, quiet)
    return CallbackSet(
        termination_forwards(typeof(pnsystem), vₑ, quiet),
        dtmin_terminator(typeof(pnsystem), eltype(pnsystem), quiet),
        decreasing_v_terminator(typeof(pnsystem), quiet),
        nonfinite_terminator(),
    )
end

function default_termination_criteria_backwards(pnsystem, v₁, quiet)
    return CallbackSet(
        termination_backwards(typeof(pnsystem), v₁, quiet),
        dtmin_terminator(pnsystem, eltype(pnsystem), quiet),
        nonfinite_terminator(),
    )
end

function default_reltol(pnsystem)
    T = eltype(pnsystem)
    return eps(T)^(11//16)
end

function default_abstol(pnsystem)
    T = eltype(pnsystem)
    return T[
        T[eps(T(M₁(pnsystem) + M₂(pnsystem)))^(11//16) for _ ∈ 1:2]
        T[eps(T)^(11//16) for _ ∈ 3:length(pnsystem.state)]
    ]
end

"""
    orbital_evolution(pnsystem; kwargs...)
    orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; kwargs...)

Integrate the orbital dynamics of an inspiraling non-eccentric compact binary.


## Required arguments

The first argument to this function may be a single `PNSystem` that encodes these required
arguments (as well as `Rᵢ`, `Λ₁`, and `Λ₂` among the keyword arguments), or the following
may be given explicitly:

  * `M₁`: Initial mass of object 1
  * `M₂`: Initial mass of object 2
  * `χ⃗₁`: Initial dimensionless spin of object 1, `S⃗₁/M₁²`
  * `χ⃗₂`: Initial dimensionless spin of object 2, `S⃗₂/M₂²`
  * `Ωᵢ`: Initial orbital angular frequency

(Note that the explicit inputs require `Ωᵢ`, whereas `PNSystem`s require `vᵢ` as input.)

These parameters all describe the "initial" conditions.  See below for an explanation of the
different meanings of "initial" and "first" in this context.  Note that the masses change in
time as a result of tidal heating — though the changes are quite small throughout most of
the inspiral.  The spins change direction due to precession, but also change in magnitude
due to tidal heating.  Therefore, the values passed here are only precisely as given
*precisely at* the moment of the initial data corresponding to the frequency `Ωᵢ`.


## Keyword arguments

Note that several of these keywords are given as Unicode but can also be given as the ASCII
string noted.  For example, `Λ₁` may be input as `Lambda1` equivalently; the default values
are the same, regardless.

  * `Λ₁=0` or `Lambda1`: Tidal-coupling parameter of object 1.
  * `Λ₂=0` or `Lambda2`: Tidal-coupling parameter of object 2.
  * `Ω₁=Ωᵢ` or `Omega_1`: First angular frequency in output data.  This may be less than
    `Ωᵢ`, in which case we integrate backwards to this point, and combine the backwards and
    forwards solutions into one seamless output.  (See next section.)
  * `Ωₑ=Ω(v=1,M=M₁+M₂)` or `Omega_e`: Final angular frequency at which to stop ODE
    integration.  Note that integration may stop before the system reaches this frequency,
    if we detect that PN has broken down irretrievably — for example, if one of the masses
    is no longer strictly positive, if a spin is super-extremal, or the PN velocity
    parameter `v` is decreasing, or is no longer in the range `(0,1)`.  Warnings will
    usually only be issued if `v < 0.35`, but if `quiet=true` informational messages will be
    issued.
  * `Rᵢ=Rotor(1)` or `R_i`: Initial orientation of binary.
  * `approximant="TaylorT1"`: Method of evaluating the right-hand side of the evolution
    equations.  Other possibilities are [`"TaylorT4"`](@ref TaylorT4!) and
    [`"TaylorT5"`](@ref TaylorT5!).  See the documentation of [`TaylorT1!`](@ref) for more
    details.
  * `PNOrder=typemax(Int)`: Order to which to retain powers of ``v^2`` in PN expansions.
    The default is to include all available terms in each PN expression.
  * `check_up_down_instability=true`: Warn if the "up-down instability" (see below) is
    likely to affect this system.
  * `time_stepper=Vern9()`: Choice of solver in OrdinaryDiffEq to integrate ODE.
  * `abstol=eps(T)^(11//16)`: Absolute tolerance of ODE solver, where `T` is the common type
    to which all the positional arguments are promoted.  This is the tolerance on local
    error estimates, not necessarily the global error.  Note that `11//16` is just chosen to
    suggest that we will have roughly 11 digits of accuracy (locally) for `Float64`
    computations, and a similar accuracy for other float types *relative to* that type's
    epsilon.
  * `reltol=eps(T)^(11//16)`: Relative tolerance of ODE solver.  (As above.)
  * `termination_criteria_forwards=nothing`: Callbacks to use for forwards-in-time
    evolution.  See below for discussion of the default value.
  * `termination_criteria_backwards=nothing`: Callbacks to use for backwards-in-time
    evolution.  See below for discussion of the default value.
  * `force_dtmin=true`: If `dt` decreases below the integrator's own minimum, and this is
    false, the integrator will immediately raise an error, before the termination criteria
    have the chance to exit gracefully.  Note that a true value here is critical if the
    `dtmin_terminator` callback is to have any effect.
  * `quiet=true`: If set to `false`, informational messages about successful terminations of
    the ODE integrations (which occur when the target ``v`` is reached in either direction)
    will be provided.  Warnings will still be issued when terminating for other reasons; if
    you wish to silence them too, you should do something like
    ```julia
    using Logging
    with_logger(SimpleLogger(Logging.Error)) do
        <your code goes here>
    end
    ```
  * `saves_per_orbit=0`: If greater than 0, the output will be interpolated so that there
    are `saves_per_orbit` time steps in the output for each orbit.  Note that this conflicts
    with the `saveat` option noted below.

All remaining keyword arguments are passed to the [`solve`
function](https://github.com/SciML/DiffEqBase.jl/blob/8e6173029c630f6908252f3fc28a69c1f0eab456/src/solve.jl#L393)
of `DiffEqBase`.  See that function's documentation for details, including useful keyword
arguments.  The most likely important one is

  * `saveat`: Denotes specific times to save the solution at, during the solving phase —
    either a time step or a vector of specific times.

In particular, if you want the solution to be output at uniform time steps `δt`, you want to
pass something like `saveat=δt`; you *don't want* the `solve` keyword `dt`, which is just
the initial suggestion for adaptive systems.  It is not permitted to pass this option *and*
the `saves_per_orbit` option.

Also note that `callback` can be used, and is combined with the callbacks generated by the
`termination_criteria_*` arguments above.  That is, you can use the default ones *and* your
own by passing arguments to `callback`.  See [the
documentation](https://diffeq.sciml.ai/dev/features/callback_functions/) for more details,
but note that if you want to make your own callbacks, you will need to add `OrdinaryDiffEq`
to your project — or possibly even `DifferentialEquations` for some of the fancier built-in
callbacks.


## ODE system

The evolved variables, in order, are

  * `M₁`: Mass of black hole 1
  * `M₂`: Mass of black hole 2
  * `χ⃗₁ˣ`: ``x`` component of dimensionless spin of black hole 1
  * `χ⃗₁ʸ`: ``y`` component...
  * `χ⃗₁ᶻ`: ``z`` component...
  * `χ⃗₂ˣ`: ``x`` component of dimensionless spin of black hole 2
  * `χ⃗₂ʸ`: ``y`` component...
  * `χ⃗₂ᶻ`: ``z`` component...
  * `Rʷ`: Scalar component of frame rotor
  * `Rˣ`: ``x`` component...
  * `Rʸ`: ``y`` component...
  * `Rᶻ`: ``z`` component...
  * `v`: PN "velocity" parameter related to the total mass ``M`` and orbital angular
    frequency ``Ω`` by ``v = (M Ω)^{1/3}``
  * `Φ`: Orbital phase given by integrating ``Ω``

The masses and spin magnitudes evolve according to [`tidal_heating`](@ref).  The spin
directions evolve according to [`Ω⃗ᵪ₁`](@ref) and [`Ω⃗ᵪ₂`](@ref).  The frame precesses with
angular velocity [`Ω⃗ₚ`](@ref), while also rotating with angular frequency `Ω` about the
[Newtonian orbital angular velocity direction](@ref ℓ̂).  The frame rotor ``R`` is given by
integrating the sum of these angular velocities as described in [Boyle
(2016)](https://arxiv.org/abs/1604.08139).  And finally, the PN parameter ``v`` evolves
according to something like
```math
\\dot{v} = - \\frac{\\mathcal{F} + \\dot{M}_1 + \\dot{M}_2} {\\mathcal{E}'}
```
where [`𝓕`](@ref) is the flux of gravitational-wave energy out of the system,
``\\dot{M}_1`` and ``\\dot{M}_2`` are due to tidal coupling as computed by
[`tidal_heating`](@ref), and [`𝓔′`](@ref) is the derivative of the binding energy with
respect to ``v``.  For `"TaylorT1"`, the right-hand side of this equation is evaluated as
given; for `"TaylorT4"`, the right-hand side is first expanded as a Taylor series in ``v``
and then truncated at some desired order; for `"TaylorT5"`, the *inverse* of the right-hand
side is expanded as a Taylor series in ``v``, truncated at some desired order, and then
inverted to obtain an expression in terms of ``v``.


## Returned solution

The returned quantity is an [`ODESolution`](https://diffeq.sciml.ai/dev/basics/solution/)
object, which has various features for extracting and interpolating the data.  We'll call
this object `sol`.

!!! note

    The solution comes with data at the time points the ODE integrator happened to
    step to.  However, it *also* comes with dense output (unless you manually turn it
    off when calling `orbital_evolution`).  This means that you can interpolate the
    solution to any other set of time points you want simply by calling it as
    `sol(t)` for some vector of time points `t`.  The quantity returned by that will
    have all the features described below, much like the original solution.  Note
    that if you only want some of the data, you can provide the optional keyword
    argument `idxs` to specify which of the elements described below you want to
    interpolate.  For example, if you only want to interpolate the values of `M₁` and
    `M₂`, you can use `sol(t, idxs=[1,2])`.

The field `sol.t` is the set of time points at which the solution is given.  To access the
`i`th variable at time step `j`, use `sol[i, j]`.[^1] You can also use colons.  For example,
`sol[:, j]` is a vector of all the data at time step `j`, and `sol[i, :]` is a vector of the
`i`th variable at all times.

[^1]: Here, the `i`th variable just refers to which number it has in the list of evolved
      variables in the ODE system, as described under "ODE system".

For convenience, you can also access the individual variables with their symbols.  For
example, `sol[:v]` returns a vector of the PN velocity parameter at each time step.  Note
the colon in `:v`, which is [Julia's notation for a
`Symbol`](https://docs.julialang.org/en/v1/base/base/#Core.Symbol).

## Initial frequency vs. first frequency vs. end frequency

Note the distinction between `Ωᵢ` (with subscript `i`) and `Ω₁` (with subscript `1`).  The
first, `Ωᵢ`, represents the angular frequency of the *initial condition* from which the ODE
integrator will begin; the second, `Ω₁`, represents the target angular frequency of the
first element of the output data.  That is, the ODE integration will run forwards in time
from `Ωᵢ` to the merger, and then — if `Ωᵢ>Ω₁` — come back to `Ωᵢ` and run backwards in time
to `Ω₁`.  The output data will stitch these two together to be one continuous
(forwards-in-time) data series.

For example, if you are trying to match to a numerical relativity (NR) simulation, you can
read the masses and spins off of the NR data when the system is orbiting at angular
frequency `Ωᵢ`.  Integrating the post-Newtonian (PN) solution forwards in time from this
point will allow you to compare the PN and NR waveforms.  However, you may want to know what
the waveform was at *earlier* times than are present in the NR data.  For this, you also
have to integrate backwards in time.  We parameterize the point to which you integrate
backwards with `Ω₁`.  In either case, element `1` of the output solution will have frequency
`Ω₁` — though by default it is equal to `Ωᵢ`.

Similarly, the optional argument `Ωₑ=1` is the frequency of the `end` element of the
solution — that is Julia's notation for the last element.  Note that this is automatically
reduced if necessary so that the corresponding PN parameter ``v`` is no greater than 1,
which may be the case whenever the total mass is greater than 1.


## Up-down instability

Be aware that the [up-down instability](http://arxiv.org/abs/1506.09116) (where the more
massive black hole has spin aligned with the orbital angular velocity, and the less massive
has spin anti-aligned) can cause systems with nearly zero precession at the initial time to
evolve into a highly precessing system either at earlier or later times.  This is a real
physical result, rather than a numerical issue.  If you want to simulate a truly
non-precessing system, you should explicitly set the in-place components of spin to
precisely 0.  By default, we check for this condition, and will issue a warning if it is
likely to be encountered for systems with low initial precession.  The function used to
compute the unstable region is [`up_down_instability`](@ref).


## Time-stepper algorithms

`Tsit5()` is a good default choice for time stepper when using `Float64` with medium-low
tolerance.  If stiffness seems to be impacting the results, `AutoTsit5(Rosenbrock23())` will
automatically switch when stiffness occurs.  For tighter tolerances, especially when using
`Double64`s, `Vern9()` or `AutoVern9(Rodas5P())` are good choices.  For very loose
tolerances, as when using `Float32`s, it might be better to use `OwrenZen3()`.


## Termination criteria

The termination criteria are vital to efficiency of the integration and correctness of the
solution.  The default values for forwards- and backwards-in-time evolution, respectively,
are
```julia
CallbackSet(
    termination_forwards(v(Ω=Ωₑ, M=M₁+M₂)),
    dtmin_terminator(T),
    decreasing_v_terminator(),
    nonfinite_terminator()
)
```
and
```julia
CallbackSet(
    termination_backwards(v(Ω=Ω₁, M=M₁+M₂)),
    dtmin_terminator(T),
    nonfinite_terminator()
)
```
where `T` is the common float type of the input arguments.  If any additional termination
criteria are needed, they could be added as additional elements of the `CallbackSet`s.  See
the [callback documentation](https://diffeq.sciml.ai/stable/features/callback_functions/)
for details.
"""
Base.@constprop :aggressive function orbital_evolution(
    M₁,
    M₂,
    χ⃗₁,
    χ⃗₂,
    Ωᵢ;
    Lambda1=0,
    Lambda2=0,
    Omega_1=0,
    Omega_e=Ω(; v=1, M=M₁ + M₂),
    R_i=Rotor(true),
    Λ₁=Lambda1,
    Λ₂=Lambda2,
    Ω₁=Omega_1,
    Ωₑ=Omega_e,
    Rᵢ=R_i,
    approximant="TaylorT1",
    PNOrder=max_pn_order,
    check_up_down_instability=true,
    time_stepper=Vern9(),
    reltol=nothing,
    abstol=nothing,
    termination_criteria_forwards=nothing,
    termination_criteria_backwards=nothing,
    quiet=true,
    force_dtmin=true,
    saves_per_orbit=false,
    solve_kwargs...,
)
    # Sanity checks for the inputs

    if M₁ ≤ 0 || M₂ ≤ 0
        error("Unphysical masses: M₁=$M₁, M₂=$M₂.")
    end

    χ⃗₁, χ⃗₂ = QuatVec(χ⃗₁), QuatVec(χ⃗₂)
    if abs2vec(χ⃗₁) > 1 || abs2vec(χ⃗₂) > 1
        error(
            "Unphysical spins: |χ⃗₁|=$(abs2vec(χ⃗₁)), |χ⃗₂|=$(abs2vec(χ⃗₂)).\n" *
            "These are dimensionless spins, which should be less than 1.\n" *
            "Perhaps you forgot to divide by M₁² or M₂², respectively.",
        )
    end

    Rᵢ = Rotor(Rᵢ)

    vᵢ = v(; Ω=Ωᵢ, M=M₁ + M₂)
    if vᵢ ≥ 1
        error(
            "The input Ωᵢ=$Ωᵢ is too large; with these masses, it corresponds to\n" *
            "vᵢ=$vᵢ, which is beyond the reach of post-Newtonian methods.",
        )
    end

    if !iszero(Λ₁) && iszero(Λ₂)
        error(
            "By convention, the NS in a BHNS binary must be the second body,\n" *
            "meaning that Λ₁ should be zero, and only Λ₂ should be nonzero.\n" *
            "You may want to swap the masses, spins, and Λ parameters.\n" *
            "Alternatively, both can be nonzero, resulting in an NSNS binary.",
        )
    end

    if Ω₁ > Ωᵢ
        error(
            "Initial frequency Ωᵢ=$Ωᵢ should be greater than " *
            "or equal to first frequency Ω₁=$Ω₁.",
        )
    end

    if Ωᵢ > Ωₑ
        error(
            "Initial frequency Ωᵢ=$Ωᵢ should be less than " *
            "or equal to ending frequency Ωₑ=$Ωₑ.",
        )
    end

    if saves_per_orbit != false && "saveat" ∈ keys(solve_kwargs)
        error(
            "It doesn't make sense to pass the `saves_per_orbit` argument *and* the " *
            "`saveat` argument; only one may be passed.",
        )
    end

    v₁ = v(; Ω=Ω₁, M=M₁ + M₂)
    vₑ = min(v(; Ω=Ωₑ, M=M₁ + M₂), 1)
    Φ = 0

    # Initial conditions for the ODE integration
    pnsystem = let R = Rᵢ, v = vᵢ
        if !iszero(Λ₁) && !iszero(Λ₂)
            NSNS(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Λ₁, Λ₂, Φ, PNOrder)
        elseif !iszero(Λ₂)
            BHNS(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Λ₂, Φ, PNOrder)
        else
            BBH(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        end
    end

    RHS! = if approximant == "TaylorT1"
        TaylorT1RHS(typeof(pnsystem))
    elseif approximant == "TaylorT4"
        TaylorT4RHS(typeof(pnsystem))
    elseif approximant == "TaylorT5"
        TaylorT5RHS(typeof(pnsystem))
    else
        error("Approximant `$approximant` is not currently supported")
    end

    if isnothing(termination_criteria_forwards)
        termination_criteria_forwards = default_termination_criteria_forwards(
            pnsystem, vₑ, quiet
        )
    end

    if isnothing(termination_criteria_backwards) && v₁ < v(pnsystem)
        termination_criteria_backwards = default_termination_criteria_backwards(
            pnsystem, v₁, quiet
        )
    end

    # The choice of 11//16 here is just an easy way to get an idea that for Float64 this
    # will give us around 11 digits of accuracy, and a similar fraction of the precision for
    # other types.
    T = eltype(pnsystem)
    if isnothing(reltol)
        reltol = default_reltol(pnsystem)
    end
    if isnothing(abstol)
        abstol = default_abstol(pnsystem)
    end

    return orbital_evolution(
        pnsystem;
        RHS!,
        v₁,
        vₑ,
        check_up_down_instability,
        quiet,
        termination_criteria_forwards,
        termination_criteria_backwards,
        time_stepper,
        reltol,
        abstol,
        force_dtmin,
        saves_per_orbit,
        solve_kwargs...,
    )
end

Base.@constprop :aggressive function orbital_evolution(
    pnsystemᵢ;
    (RHS!)=(TaylorT1RHS(typeof(pnsystemᵢ))),
    v₁=zero(pnsystemᵢ),
    vₑ=one(pnsystemᵢ),
    check_up_down_instability=true,
    quiet=true,
    termination_criteria_forwards=default_termination_criteria_forwards(
        pnsystemᵢ, vₑ, quiet
    ),
    termination_criteria_backwards=default_termination_criteria_backwards(
        pnsystemᵢ, v₁, quiet
    ),
    time_stepper=Vern9(),
    reltol=default_reltol(pnsystemᵢ),
    abstol=default_abstol(pnsystemᵢ),
    force_dtmin=true,
    saves_per_orbit=zero(eltype(pnsystemᵢ)),
    solve_kwargs...,
)
    pnsystem = deepcopy(pnsystemᵢ)

    if check_up_down_instability
        up_down_instability_warn(pnsystem, v₁, vₑ)
    end

    # Log an error if the initial parameters return a NaN on the right-hand side
    let
        uᵢ = copy(pnsystem.state)
        u̇ = similar(uᵢ)
        tᵢ = zero(eltype(pnsystem))
        RHS!(u̇, uᵢ, pnsystem, tᵢ)
        if any(isnan, u̇) || any(isnan, uᵢ)
            # COV_EXCL_START
            @error "Found a NaN with initial parameters:" value.(uᵢ) value.(u̇) pnsystem
            error("Found NaN")
            # COV_EXCL_STOP
        end
    end

    # Note: This estimate for the time span over which to integrate may be very bad,
    # especially close to merger.  An underestimate would lead to an inspiral ending too
    # soon, but an overestimate can lead to integration continuing very slowly in a regime
    # where PN has broken down.
    tₑ = 4estimated_time_to_merger(pnsystem)
    tᵢ = zero(tₑ)
    if "saveat" ∈ keys(solve_kwargs) && solve_kwargs["saveat"] isa AbstractVector
        tₑ = max(tᵢ, min(tₑ, solve_kwargs["saveat"][end]))
    end

    problem_forwards = ODEProblem(
        RHS!, pnsystem.state, (tᵢ, tₑ), pnsystem; callback=termination_criteria_forwards
    )

    solution_forwards = solve(
        problem_forwards, time_stepper; reltol, abstol, force_dtmin, solve_kwargs...
    )

    solution = if v₁ > 0
        # Note: Here again, we don't want to overestimate the time span by too much, but we
        # also don't want to underestimate and get a shortened waveform.  This should be a
        # better estimate, though, because it's dealing with lower speeds, at which PN
        # approximation should be more accurate.
        pnsystem.state[:] .= pnsystemᵢ.state
        # NOTE this only makes sense for PN systems where v is a state
        # variable. That's not appropriate for eccentric PN, where it would be
        # better to use x as a state variable (or n)
        pnsystem[:v] = v₁
        t₁ =
            -4 * (estimated_time_to_merger(pnsystem) - estimated_time_to_merger(pnsystemᵢ))
        if "saveat" ∈ keys(solve_kwargs) && solve_kwargs["saveat"] isa AbstractVector
            t₁ = min(tᵢ, max(t₁, solve_kwargs["saveat"][begin]))
        end

        # Reset state to initial conditions
        pnsystem.state[:] .= pnsystemᵢ.state

        problem_backwards = remake(
            problem_forwards; tspan=(tᵢ, t₁), callback=termination_criteria_backwards
        )

        solution_backwards = solve(
            problem_backwards,
            time_stepper;
            reltol,
            abstol,
            force_dtmin,
            solve_kwargs...,
        )

        combine_solutions(solution_backwards, solution_forwards)
    else
        solution_forwards
    end

    if saves_per_orbit > 0
        solution = uniform_in_phase(solution, saves_per_orbit)
    end

    return solution
end
