using Base: promote_typeof

"""
    value(x)

Return `x` or the value wrapped by the `Dual` number `x`

"""
value(x) = hasproperty(x, :value) ? getproperty(x, :value) : x


"""
    termination_forwards(vₑ, [quiet])

Construct termination criteria of solving PN evolution forwards in time

These criteria include checking that the masses are positive and the
dimensionless spins are less than 1, as well as ensuring that the evolution
will terminate at `vₑ`.

The optional `quiet` argument will silence informational messages about
reaching the target value of `vₑ` if set to `true`, but warnings will still be
issued when terminating for other reasons.  If you want to quiet warnings also,
you can do something like this:
```julia
using Logging
with_logger(SimpleLogger(Logging.Error)) do
    <your code goes here>
end
```

"""
function termination_forwards(vₑ, quiet=false)
    # Triggers the `continuous_terminator!` whenever one of these conditions crosses 0.
    # More precisely, the integrator performs a root find to finish precisely
    # when one of these conditions crosses 0.
    function conditions(out,u,t,integrator)
        out[1] = u[1]  # Terminate if M₁≤0
        out[2] = u[2]  # Terminate if M₂≤0
        out[3] = 1 - abs2vec(QuatVec{typeof(vₑ)}(u[3:5]...))  # Terminate if χ₁>1
        out[4] = 1 - abs2vec(QuatVec{typeof(vₑ)}(u[6:8]...))  # Terminate if χ₂>1
        out[5] = vₑ - u[13]  # Terminate at v = vₑ
    end
    function terminator!(integrator, event_index)
        if event_index == 1
            @warn "Terminating forwards evolution because M₁ has become non-positive.  This is unusual."
        elseif event_index == 2
            @warn "Terminating forwards evolution because M₂ has become non-positive.  This is unusual."
        elseif event_index == 3
            @warn "Terminating forwards evolution because χ₁>1.  Suggests early breakdown of PN."
        elseif event_index == 4
            @warn "Terminating forwards evolution because χ₂>1.  Suggests early breakdown of PN."
        elseif event_index == 5
            quiet || @info (
                "Terminating forwards evolution because the PN parameter 𝑣 "
                * "has reached 𝑣ₑ=$(value(vₑ)).  This is ideal."
            )
        end
        terminate!(integrator)
    end
    VectorContinuousCallback(
        conditions,
        terminator!,
        5;  # We have 5 criteria above
        save_positions=(false,false)  # Only save before the termination, not after
    )
end


"""
    termination_backwards(v₁, [quiet])

Construct termination criteria of solving PN evolution backwards in time

These criteria include checking that the masses are positive and the
dimensionless spins are less than 1, as well as ensuring that the evolution
will terminate at `v₁`.

The optional `quiet` argument will silence informational messages about
reaching the target value of `v₁` if set to `true`, but warnings will still be
issued when terminating for other reasons.  If you want to quiet warnings also,
you can do something like this:
```julia
using Logging
with_logger(SimpleLogger(Logging.Error)) do
    <your code goes here>
end
```

"""
function termination_backwards(v₁, quiet=false)
    function terminators_backwards(out,u,t,integrator)
        out[1] = u[1]  # Terminate if M₁≤0
        out[2] = u[2]  # Terminate if M₂≤0
        out[3] = 1 - abs2vec(QuatVec{typeof(v₁)}(u[3:5]...))  # Terminate if χ₁>1
        out[4] = 1 - abs2vec(QuatVec{typeof(v₁)}(u[6:8]...))  # Terminate if χ₂>1
        out[5] = v₁ - u[13]  # Terminate at v = v₁
    end
    function terminator_backwards!(integrator, event_index)
        if event_index == 1
            @warn "Terminating backwards evolution because M₁ has become non-positive.  Suggests problem with PN."
        elseif event_index == 2
            @warn "Terminating backwards evolution because M₂ has become non-positive.  Suggests problem with PN."
        elseif event_index == 3
            @warn "Terminating backwards evolution because χ₁>1.  Suggests problem with PN."
        elseif event_index == 4
            @warn "Terminating backwards evolution because χ₂>1.  Suggests problem with PN."
        elseif event_index == 5
            quiet || @info (
                "Terminating backwards evolution because the PN parameter 𝑣 "
                * "has reached 𝑣₁=$(value(v₁)).  This is ideal."
            )
        end
        terminate!(integrator)
    end
    VectorContinuousCallback(
        terminators_backwards,
        terminator_backwards!,
        5;  # We have 5 criteria above
        save_positions=(false,false)  # Only save before the termination, not after
    )
end


"""
    dtmin_terminator(T)

Construct termination criterion to terminate when `dt` drops below `√eps(T)`.

Pass `force_dtmin=true` to `solve` when using this callback.  Otherwise, the
time-step size may decrease too much *within* a single time step, so that the
integrator itself will quit before reaching this callback, leading to a less
graceful exit.

"""
function dtmin_terminator(T)
    # Triggers the `discrete_terminator!` whenever this condition is true after
    # an integration step
    ϵ = √eps(T)
    function discrete_condition(u,t,integrator)
        abs(integrator.dt) < ϵ
    end
    function discrete_terminator!(integrator)
        @warn "Terminating forwards evolution because |dt=$(integrator.dt)| < ϵ=$(ϵ)"
        terminate!(integrator)
    end
    DiscreteCallback(
        discrete_condition,
        discrete_terminator!;
        save_positions=(false,false)
    )
end


"""
    nonfinite_terminator()

Construct termination criterion to terminate when any NaN or Inf is found in
the data after an integration step.

"""
function nonfinite_terminator()
    # Triggers the `discrete_terminator!` whenever this condition is true after
    # an integration step
    function discrete_condition(u,t,integrator)
        # any(isnan, u) || isnan(t) || isnan(integrator.dt)
        !(all(isfinite, u) && isfinite(t) && isfinite(integrator.dt))
    end
    function discrete_terminator!(integrator)
        @warn "Terminating forwards evolution because a non-finite number was found"
        terminate!(integrator)
    end
    DiscreteCallback(
        discrete_condition,
        discrete_terminator!;
        save_positions=(false,false)
    )
end


"""
    inspiral(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; kwargs...)

Integrate the orbital dynamics of an inspiraling non-eccentric compact binary.


## Keyword arguments

  * `Ω₁=Ωᵢ`: First angular velocity in output data (see next section).
  * `Ωₑ=1`: Final angular velocity at which to stop ODE integration.
  * `Rᵢ=Rotor(true)`: Initial orientation of binary.
  * `PNSys=TaylorT1`: Currently the only possibility.
  * `PNOrder=7//2`: Not actually used currently.
  * `check_up_down_instability=true`: Warn if the "up-down instability" (see
    below) is likely to affect this system.
  * `time_stepper=AutoVern9(Rodas5())`: Choice of solver in OrdinaryDiffEq to
    integrate ODE.
  * `abstol=eps(T)^(11//16)`: Absolute tolerance of ODE solver, where `T` is
    the common type to which all the positional arguments are promoted.  This
    is the tolerance on local error estimates, not necessarily the global
    error.
  * `reltol=eps(T)^(11//16)`: Relative tolerance of ODE solver.  (As above.)
  * `termination_criteria_forwards=nothing`: Callbacks to use for
    forwards-in-time evolution.  See below for discussion of the default value.
  * `termination_criteria_backwards=nothing`: Callbacks to use for
    backwards-in-time evolution.  See below for discussion of the default value.
  * `force_dtmin=true`: If `dt` decreases below the integrator's own minimum,
    and this is false, the integrator will immediately raise an error, before
    the termination criteria have the chance to exit gracefully.  Note that a
    true value here is critical if the `dtmin_terminator` callback is to have
    any effect.
  * `quiet=false`: If set to `true`, informational messages about successful
    terminations of the ODE integrations (which occur when the target ``v`` is
    reached in either direction) will be silenced.  Warnings will still be
    issued when terminating for other reasons; if you wish to silence them too,
    you should do something like
    ```julia
    using Logging
    with_logger(SimpleLogger(Logging.Error)) do
        <your code goes here>
    end
    ```
  * `integrate_orbital_phase=false`: If set to `true`, integrate the orbital
    phase ``Φ`` along with the rest of the system.  Note that this may slow the
    system down because the absolute value of ``Φ`` may grow to very large
    values, so that the `abstol` will strain to keep its evolution far more
    accurate than is really needed.  If this is a problem, you can loosen
    `abstol` and/or pass vectors of separate tolerances for each variable in
    the ODE system (see below).

All remaining keyword arguments are passed to the [`solve`
function](https://github.com/SciML/DiffEqBase.jl/blob/8e6173029c630f6908252f3fc28a69c1f0eab456/src/solve.jl#L393)
of `DiffEqBase`.  See that function's documentation for details, including
useful keyword arguments.  The most likely important ones are

  * `saveat`: Denotes specific times to save the solution at, during the
    solving phase.
  * `dt`: Sets the *initial* stepsize. Defaults to an automatic choice if the
    method is adaptive.
  * `dtmax`: Maximum dt for adaptive timestepping.
  * `dtmin`: Minimum dt for adaptive timestepping.

Note that if you want the solution to be output with specifically spaced time
steps, you *don't* want `dt`, which is just the initial suggestion for adaptive
systems; you want to set `saveat` to the desired spacing.  [The `saveat`
argument could be a vector of specific times at which to save, but because we
don't know when the PN evolution ends, this probably isn't useful.]

Also note that `callback` can be used, and is combined with the callbacks
generated by the `termination_criteria_*` arguments above.  See [the
documentation](https://diffeq.sciml.ai/dev/features/callback_functions/) for
more details, but note that if you want to make your own callbacks, you will
need to add `OrdinaryDiffEq` to your project — or possibly even
`DifferentialEquations` for some of the fancier built-in callbacks.


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
  * `v`: PN "velocity" parameter related to the total mass ``M`` and orbital
    angular velocity ``Ω`` by ``v = (M Ω)^{1/3}``
  * `Φ`: Orbital phase given by integrating ``Ω`` (optional; only appears if
    `integrate_orbital_phase` is `true`)

The masses and spin magnitudes evolve according to [`tidal_heating`](@ref).
The spin directions evolve according to [`Ω⃗ᵪ₁`](@ref) and [`Ω⃗ᵪ₂`](@ref).  The
frame rotor ``R`` is given by integrating the angular velocity as described in
[Boyle (2016)](https://arxiv.org/abs/1604.08139), while the angular velocity
itself is given by [`Ω⃗ₚ`](@ref).  And finally, the PN parameter ``v`` evolves
according to something like
```math
\\dot{v} = - \\frac{\\mathcal{F} + \\dot{M}_1 + \\dot{M}_2} {\\mathcal{E}'}
```
where [`𝓕`](@ref) is the flux of gravitational-wave energy out of the system
and [`𝓔′`](@ref) is the derivative of the binding energy with respect to ``v``.
For `"TaylorT1"`, the right-hand side of this equation is evaluated as given;
for `"TaylorT4"`, the right-hand side is first expanded as a Taylor series in
``v`` and then truncated at some desired order; for `"TaylorT5"`, the *inverse*
of the right-hand side is expanded as a Taylor series in ``v``, truncated at
some desired order, and then inverted to obtain an expression in terms of
``v``.


## Returned solution

The returned quantity is an
[`ODESolution`](https://diffeq.sciml.ai/dev/basics/solution/) object, which has
various features for extracting and interpolating the data.  We'll call this
object `sol`.

!!! note
    
    The solution comes with data at the time points the ODE integrator happened
    to step to.  However, it *also* comes with dense output (unless you
    manually turn it off when calling `inspiral`).  This means that you can
    interpolate the solution to any other set of time points you want simply by
    calling it as `sol(t)` for some vector of time points `t`.  The quantity
    returned by that will have the following features, just like the original
    solution.  Note that if you only want some of the data you can provide the
    optional keyword argument `idxs` to specify which of the elements described
    below you want to interpolate.  For example, if you only want to
    interpolate the values of `M₁` and `M₂`, you can use `sol(t, idxs=[1,2])`.

The field `sol.t` is the set of time points at which the solution is given.  To
access the `i`th variable at time step `j`, use `sol[i, j]`.[^1] You can also
use colons.  For example, `sol[:, j]` is a vector of all the data at time step
`j`, and `sol[i, :]` is a vector of the `i`th variable at all times.

[^1]: Here, the `i`th variable just refers to which number it has in the list
      of evolved variables in the ODE system, as described under "ODE system".


## Initial frequency vs. first frequency vs. end frequency

Note the distinction between `Ωᵢ` (with subscript `i`) and `Ω₁` (with subscript
`1`).  The first, `Ωᵢ`, represents the angular velocity of the *initial
condition* from which the ODE integrator will begin; the second, `Ω₁`,
represents the target angular velocity of the first element of the output data.
That is, the ODE integration will run forwards in time from `Ωᵢ` to the merger,
and then come back to `Ωᵢ` and run backwards in time to `Ω₁`.  The output data
will stitch these two together to be one continuous (forwards-in-time) data
series.

For example, if you are trying to match to a numerical relativity (NR)
simulation, you can read the masses and spins off of the NR data when the
system is orbiting at angular velocity `Ωᵢ`.  Integrating the post-Newtonian
(PN) solution forwards in time from this point will allow you to compare the PN
and NR waveforms.  However, you may want to know what the waveform was at
*earlier* times than are present in the NR data.  For this, you also have to
integrate backwards in time.  We parametrise the point to which you integrate
backwards with `Ω₁`.  In either case, element `1` of the output solution will
have frequency `Ω₁` — though by default it is equal to `Ωᵢ`.

Similarly, the optional argument `Ωₑ=1` is the frequency of the `end` element
of the solution — that is Julia's notation for the last element.  Note that
this is automatically reduced if necessary so that the corresponding PN
parameter ``v`` is no greater than 1, which may be the case whenever the total
mass is greater than 1.


## Up-down instability

Be aware that the [up-down instability](http://arxiv.org/abs/1506.09116) (where
the more massive black hole has spin aligned with the orbital angular velocity,
and the less massive has spin anti-aligned) can cause systems with nearly zero
precession at the initial time to evolve into a highly precessing system either
at earlier or later times.  This is a real physical result, rather than a
numerical issue.  If you want to simulate a truly non-precessing system, you
should explicitly set the in-place components of spin to precisely 0.  By
default, we check for this condition, and will issue a warning if it is likely
to be encountered for systems with low initial precession.  The function used
to compute the unstable region is [`up_down_instability`](@ref).


## Time-stepper algorithms

`Tsit5()` is a good default choice for time stepper when using `Float64` with
medium-low tolerance.  If stiffness seems to be impacting the results,
`AutoTsit5(Rosenbrock23())` will automatically switch when stiffness occurs.
For tighter tolerances, especially when using `Double64`s, `Vern9()` or
`AutoVern9(Rodas5())` are good choices.  For very loose tolerances, as when
using `Float32`s, it might be better to use `OwrenZen3()`.


## Termination criteria

The termination criteria are vital to efficiency of the integration and
correctness of the solution.  The default values for forwards- and
backwards-in-time evolution, respectively, are
```julia
CallbackSet(
    termination_forwards(v(Ω=Ωₑ, M=M₁+M₂)),
    dtmin_terminator(T),
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
where `T` is the common float type of the input arguments.  If any additional
termination criteria are needed, they could be added as additional elements of
the `CallbackSet`s.  See the [callback
documentation](https://diffeq.sciml.ai/stable/features/callback_functions/) for
details.

"""
function inspiral(
    M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ;
    Ω₁=Ωᵢ, Ωₑ=1, Rᵢ=Rotor(true),
    PNSys=TaylorT1, PNOrder=7//2,
    check_up_down_instability=true, time_stepper=AutoVern9(Rodas5()),
    reltol=nothing, abstol=nothing,
    termination_criteria_forwards=nothing,
    termination_criteria_backwards=nothing,
    force_dtmin=true, integrate_orbital_phase=false,
    quiet=false,
    solve_kwargs...
)
    if Ω₁ > Ωᵢ
        error(
            "Initial frequency Ωᵢ=$Ωᵢ should be greater than or equal to first frequency Ω₁=$Ω₁"
        )
    end

    vᵢ = v(Ω=Ωᵢ, M=M₁+M₂)
    if vᵢ ≥ 1
        error(
            "The input Ωᵢ=$Ωᵢ is too large; with these masses, it corresponds to "
            * "vᵢ=$vᵢ, which is beyond the reach of post-Newtonian methods."
        )
    end

    v₁ = v(Ω=Ω₁, M=M₁+M₂)
    vₑ = min(v(Ω=Ωₑ, M=M₁+M₂), 1)

    # Initial conditions for the ODE integration
    uᵢ = [M₁; M₂; χ⃗₁.components[2:4]; χ⃗₂.components[2:4]; Rᵢ.components; vᵢ]
    # We pack this up here, to get everything into the same type, and permit easier
    # passing to the other form of this function; we'll unpack it again there,
    # which makes sure everything has the same type and the function is type stable.

    T = eltype(uᵢ)
    if isnothing(reltol)
        reltol = eps(T)^(11//16)
    end
    if isnothing(abstol)
        abstol = eps(T)^(11//16)
    end

    inspiral(
        uᵢ, Ω₁, Ωₑ, v₁, vₑ,
        PNSys(PNOrder, T), T,
        check_up_down_instability, time_stepper,
        reltol, abstol,
        termination_criteria_forwards,
        termination_criteria_backwards,
        force_dtmin, integrate_orbital_phase,
        quiet;
        solve_kwargs...
    )
end

function inspiral(
    uᵢ, Ω₁, Ωₑ, v₁, vₑ,
    pn::PNSystem, T,
    check_up_down_instability, time_stepper,
    reltol, abstol,
    termination_criteria_forwards,
    termination_criteria_backwards,
    force_dtmin, integrate_orbital_phase,
    quiet;
    solve_kwargs...
)

    M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, vᵢ = uᵢ
    χ⃗₁ = QuatVec(χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ)
    χ⃗₂ = QuatVec(χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ)
    R = Quaternion(Rʷ, Rˣ, Rʸ, Rᶻ)

    if check_up_down_instability
        χₚₑᵣₚ = let n̂=n̂(R), λ̂=λ̂(R)
            √((χ⃗₁ ⋅ n̂)^2 + (χ⃗₁ ⋅ λ̂)^2 + (χ⃗₂ ⋅ n̂)^2 + (χ⃗₂ ⋅ λ̂)^2)
        end
        if χₚₑᵣₚ ≤ 1e-2
            (Ω₊, Ω₋) = up_down_instability(uᵢ)
            if Ω₁ < min(Ω₋, 1//2) && min(Ωₑ, 1//2) > Ω₊
                @warn (
                    "This system is likely to encounter the up-down instability in the\n"
                    * "frequency range (Ω₊, Ω₋)=$((Ω₊, Ω₋)).\n"
                    * "This is a true physical instability; not just a numerical issue.\n"
                    * "Despite the initial conditions containing very small precession,\n"
                    * "the system will likely evolve to have very large precession."
                )
            end
        end
    end

    estimated_time_to_merger = 5/(256ν(M₁, M₂) * T(vᵢ)^8) # Lowest-order PN time-to-merger
    tspan = (T(0), 4estimated_time_to_merger)
    if isnothing(termination_criteria_forwards)
        termination_criteria_forwards = CallbackSet(
            termination_forwards(vₑ, quiet),
            dtmin_terminator(T),
            nonfinite_terminator()
        )
    end
    problem_forwards = ODEProblem(
        noneccentric_RHS!,
        integrate_orbital_phase ? [uᵢ; zero(T)] : uᵢ,
        tspan, pn,
        callback=termination_criteria_forwards
    )

    # Log an error if the initial parameters return a NaN on the right-hand side
    let
        u̇ = similar(uᵢ)
        noneccentric_RHS!(u̇, uᵢ, pn, tspan[1])
        if any(isnan, u̇) ||  any(isnan, uᵢ) ||  any(isnan, tspan)
            @error "Found a NaN with initial parameters:" value.(uᵢ) value.(u̇) pn value.(tspan)
            error("Found NaN")
        end
    end

    solution_forwards = solve(
        problem_forwards, time_stepper;
        reltol=reltol, abstol=abstol,
        force_dtmin=force_dtmin,
        solve_kwargs...
    )

    if v₁ < vᵢ
        estimated_backwards_time = 5/(256ν(M₁, M₂) * T(v₁)^8) - estimated_time_to_merger
        tspan = (T(0), -3estimated_backwards_time)
        if isnothing(termination_criteria_backwards)
            termination_criteria_backwards = CallbackSet(
                termination_backwards(v₁, quiet),
                dtmin_terminator(T),
                nonfinite_terminator()
            )
        end
        problem_backwards = remake(problem_forwards; tspan=tspan, callback=termination_criteria_backwards)

        solution_backwards = solve(
            problem_backwards, time_stepper;
            reltol=reltol, abstol=abstol,
            force_dtmin=force_dtmin,
            solve_kwargs...
        )

        combine_solutions(solution_backwards, solution_forwards)
    else
        solution_forwards
    end
end


"""
    noneccentric_RHS!(u̇, u, p, t)

Compute the right-hand side for the orbital evolution of a non-eccentric binary

Here, `u` is the ODE state vector, which can be unpacked with
[`PNDynamicalVariables`](@ref).  The parameter `p` is currently unused, but
could be used to pass un-evolved parameters through.

"""
function noneccentric_RHS!(u̇, u, p, t)
    recalculate!(u̇, u, p)
    nothing
end
