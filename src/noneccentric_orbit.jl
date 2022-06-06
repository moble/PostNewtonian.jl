using Base: promote_typeof

"""
    noneccentric_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; kwargs...)

Integrate the orbital dynamics of a non-eccentric compact binary.

## Keyword arguments

  * `Ω₁=Ωᵢ`: First angular velocity in output data (see next section).
  * `Ωₑ=1`: Final angular velocity at which to stop ODE integration.
  * `Rᵢ=Rotor(true)`: Initial orientation of binary.
  * `PNSys=TaylorT1`: Not actually used currently.
  * `PNOrder=7//2`: Not actually used currently.
  * `check_up_down_instability=true`: Warn if the [Up-down instability](@ref)
    is likely to affect this system.
  * `time_stepper=DP8()`: Choice of solver in OrdinaryDiffEq to integrate ODE.
  * `abstol=eps(T)^(11//16)`: Absolute tolerance of ODE solver, where `T` is
    the common type to which all the positional arguments are promoted.  This
    is the tolerance on local error estimates, not necessarily the global
    error.
  * `reltol=eps(T)^(11//16)`: Relative tolerance of ODE solver.  (As above.)

All remaining keyword arguments are passed to the [`solve`
function](https://github.com/SciML/DiffEqBase.jl/blob/8e6173029c630f6908252f3fc28a69c1f0eab456/src/solve.jl#L393)
of `DiffEqBase`.  See that function's documentation for details, including
useful keyword arguments.  The most likely important ones are

  * `saveat`: Denotes specific times to save the solution at, during the
    solving phase.
  * `adaptive`: Turns on adaptive timestepping for appropriate methods. Default
    is true.
  * `dt`: Sets the initial stepsize. Defaults to an automatic choice if the
    method is adaptive.
  * `dtmax`: Maximum dt for adaptive timestepping.
  * `dtmin`: Minimum dt for adaptive timestepping.

Note that `callback` is already used by this function (in addition to the
`abstol` and `reltol` mentioned above), which currently makes it impossible to
modify the callbacks.  Hacking will be required to change that.

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
of the solution — that is Julia's notation for the last element.

## Up-down instability

Be aware that the [up-down instability](http://arxiv.org/abs/1506.09116) (where
the more massive black hole has spin aligned with the orbital angular velocity,
and the less massive has spin anti-aligned) can cause systems with nearly zero
precession at the initial time to evolve into a highly precessing system either
at earlier or later times.  This is a real physical result, rather than a
numerical issue.  If you want to simulate a truly non-precessing system, you
should explicitly set the in-place components of spin to precisely 0.  By
default, we check for this condition, and will issue a warning if it is likely
to be encountered for systems with low initial precession.

## Time-stepper algorithms

`Tsit5()` is a good default choice for time stepper when using `Float64`s.  If
stiffness seems to be impacting the results, `AutoTsit5(Rosenbrock23())` will
automatically switch when stiffness occurs.  For tighter tolerances, as when
using `Double64`s, `Vern9()` or `AutoVern9(Rodas5())` are good choices.  For
very loose tolerances, as when using `Float32`s, it might be better to use
`OwrenZen3()`.

"""
function noneccentric_evolution(
    M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Ω₁=Ωᵢ, Ωₑ=1, Rᵢ=Rotor(true),
    PNSys=TaylorT1, PNOrder=7//2,
    check_up_down_instability=true, time_stepper=DP8(),
    reltol=eps(promote_typeof(M₁, M₂, χ⃗₁.vec..., χ⃗₂.vec..., Ωᵢ^(1//3), Rᵢ.components...))^(11//16),
    abstol=eps(promote_typeof(M₁, M₂, χ⃗₁.vec..., χ⃗₂.vec..., Ωᵢ^(1//3), Rᵢ.components...))^(11//16),
    solve_kwargs...
)
    if Ω₁ > Ωᵢ
        error(
            "Initial frequency Ωᵢ=$Ωᵢ should be greater than or equal to first frequency Ω₁=$Ω₁"
        )
    end

    vᵢ = v(Ω=Ωᵢ, M=M₁+M₂)
    v₁ = v(Ω=Ω₁, M=M₁+M₂)
    vₑ = v(Ω=Ωₑ, M=M₁+M₂)
    uᵢ = [  # Initial conditions for the ODE integration
        M₁;
        M₂;
        χ⃗₁.vec;
        χ⃗₂.vec;
        Rᵢ.components;
        vᵢ
    ]
    T = eltype(uᵢ)
    pn = PNSys(PNOrder, T)
    unpack!(pn, uᵢ)

    if check_up_down_instability
        χₚₑᵣₚ = let n̂=n̂(pn.R),λ̂=λ̂(pn.R)
            √((pn.χ⃗₁ ⋅ n̂)^2 + (pn.χ⃗₁ ⋅ λ̂)^2 + (pn.χ⃗₂ ⋅ n̂)^2 + (pn.χ⃗₂ ⋅ λ̂)^2)
        end
        if χₚₑᵣₚ ≤ 1e-2
            (Ω₊, Ω₋) = up_down_instability(pn)
            if Ω₁ < Ω₋ < 1//4 || Ω₁ < Ω₊ < 1//4
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
    @warn "Need to complete orbital_dynamics in recalculate!"
    problem_forwards = ODEProblem(noneccentric_RHS!, uᵢ, tspan, pn)

    # Triggers the `terminator!` whenever one of these conditions crosses 0.
    # More precisely, the integrator performs a root find to finish precisely
    # when one of these conditions crosses 0.
    function terminators(out,u,t,integrator)
        out[1] = u[1]  # Terminate if M₁≤0
        out[2] = u[2]  # Terminate if M₂≤0
        out[3] = 1 - abs2vec(QuatVec{T}(u[3:5]...))  # Terminate if χ₁>1
        out[4] = 1 - abs2vec(QuatVec{T}(u[6:8]...))  # Terminate if χ₂>1
        out[5] = vₑ - u[end]  # Terminate at v = vₑ
    end
    function terminator!(integrator, event_index)
        if event_index == 1
            @info "Terminating forwards evolution because M₁ has become non-positive.  This is unusual."
        elseif event_index == 2
            @info "Terminating forwards evolution because M₂ has become non-positive.  This is unusual."
        elseif event_index == 3
            @info "Terminating forwards evolution because χ₁>1.  Suggests early breakdown of PN."
        elseif event_index == 4
            @info "Terminating forwards evolution because χ₂>1.  Suggests early breakdown of PN."
        elseif event_index == 5
            @info (
                "Terminating forwards evolution because the PN parameter 𝑣 "
                * "has reached 𝑣ₑ=$vₑ.  This is ideal."
            )
        end
        terminate!(integrator)
    end
    termination_criteria = VectorContinuousCallback(
        terminators,
        terminator!,
        5;  # We have 5 criteria above
        save_positions=(true,false)  # Only save before the termination, not after
    )

    solution_forwards = solve(
        problem_forwards, time_stepper,
        callback=termination_criteria,
        reltol=reltol, abstol=abstol,
        solve_kwargs...
    )

    if v₁ < vᵢ
        estimated_backwards_time = 5/(256ν(M₁, M₂) * T(v₁)^8) - estimated_time_to_merger
        tspan = (T(0), -3estimated_backwards_time)

        problem_backwards = remake(problem_forwards; tspan=tspan)

        function terminators_backwards(out,u,t,integrator)
            out[1] = u[1]  # Terminate if M₁≤0
            out[2] = u[2]  # Terminate if M₂≤0
            out[3] = 1 - abs2vec(QuatVec{T}(u[3:5]...))  # Terminate if χ₁>1
            out[4] = 1 - abs2vec(QuatVec{T}(u[6:8]...))  # Terminate if χ₂>1
            out[5] = v₁ - u[end]  # Terminate at v = v₁
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
                @info (
                    "Terminating backwards evolution because the PN parameter 𝑣 "
                    * "has reached 𝑣₁=$v₁.  This is ideal."
                )
            end
            terminate!(integrator)
        end
        termination_criteria_backwards = VectorContinuousCallback(
            terminators_backwards,
            terminator_backwards!,
            5;  # We have 5 criteria above
            save_positions=(true,false)  # Only save before the termination, not after
        )

        solution_backwards = solve(
            problem_backwards, time_stepper;
            callback=termination_criteria_backwards,
            reltol=reltol, abstol=abstol,
            solve_kwargs...
        )

        @warn "Failing to combine forwards and backwards!!!"
        return solution_backwards[end:-1:2], solution_forwards
    end

    solution_forwards
end


"""
    noneccentric_RHS!(u̇, u, p, t)

Compute the right-hand side for the orbital evolution of a non-eccentric binary

Here, `u` is the ODE state vector, which can be unpacked with
[`PNDynamicalVariables`](@ref).  The parameter `p` is currently unused, but
could be used to pass un-evolved parameters through.

"""
function noneccentric_RHS!(u̇, u, pn, t)
    recalculate!(u̇, u, pn)
    nothing
end
