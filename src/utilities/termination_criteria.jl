"""
    termination_forwards(vₑ, [quiet])

Construct termination criteria of solving PN evolution forwards in time

These criteria include checking that the masses are positive and the dimensionless spins are
less than 1, as well as ensuring that the evolution will terminate at `vₑ`.

The optional `quiet` argument will silence informational messages about reaching the target
value of `vₑ` if set to `true`, but warnings will still be issued when terminating for other
reasons.
"""
function termination_forwards(PNType, vₑ, quiet=false)
    # Triggers the `continuous_terminator!` whenever one of these conditions crosses 0.
    # More precisely, the integrator performs a root find to finish precisely
    # when one of these conditions crosses 0.

    # NOTE These are very specific to quasicircular; have to generalize for eccentric
    function conditions(out, state, t, integrator)
        out[1] = state[symbol_index(PNType, :M₁)] # Terminate if M₁ ≤ 0
        out[2] = state[symbol_index(PNType, :M₂)]  # Terminate if M₂ ≤ 0
        out[3] = 1 - sum(abs2, state[symbol_index(PNType, :χ⃗₁ˣ) : symbol_index(PNType, :χ⃗₁ᶻ)]) # Terminate if χ₁ > 1
        out[4] = 1 - sum(abs2, state[symbol_index(PNType, :χ⃗₂ˣ) : symbol_index(PNType, :χ⃗₂ᶻ)]) # Terminate if χ₂ > 1
        return out[5] = vₑ - state[symbol_index(PNType, :v)] # Terminate at v = vₑ
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
                "Terminating forwards evolution because the PN parameter 𝑣 " *
                "has reached 𝑣ₑ=$(value(vₑ)).  This is ideal."
            )
        end
        return terminate!(integrator)
    end
    return VectorContinuousCallback(
        conditions,
        terminator!,
        5;  # We have 5 criteria above
        save_positions=(false, false),  # Only save before the termination, not after
    )
end

"""
    termination_backwards(v₁, [quiet])

Construct termination criteria of solving PN evolution backwards in time

These criteria include checking that the masses are positive and the dimensionless spins are
less than 1, as well as ensuring that the evolution will terminate at `v₁`.

The optional `quiet` argument will silence informational messages about reaching the target
value of `v₁` if set to `true`, but warnings will still be issued when terminating for other
reasons.
"""
function termination_backwards(PNType, v₁, quiet=false)
    function terminators_backwards(out, state, t, integrator)
        out[1] = state[symbol_index(PNType, :M₁)]  # Terminate if M₁ ≤ 0
        out[2] = state[symbol_index(PNType, :M₂)]  # Terminate if M₂ ≤ 0
        out[3] = 1 - sum(abs2, state[symbol_index(PNType, :χ⃗₁ˣ) : symbol_index(PNType, :χ⃗₁ᶻ)])  # Terminate if χ₁ > 1
        out[4] = 1 - sum(abs2, state[symbol_index(PNType, :χ⃗₂ˣ) : symbol_index(PNType, :χ⃗₂ᶻ)])  # Terminate if χ₂ > 1
        return out[5] = v₁ - state[symbol_index(PNType, :v)]  # Terminate at v = v₁
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
                "Terminating backwards evolution because the PN parameter 𝑣 " *
                "has reached 𝑣₁=$(value(v₁)).  This is ideal."
            )
        end
        return terminate!(integrator)
    end
    return VectorContinuousCallback(
        terminators_backwards,
        terminator_backwards!,
        5;  # We have 5 criteria above
        save_positions=(false, false),  # Only save before the termination, not after
    )
end

"""
    dtmin_terminator(T, [quiet])

Construct termination criterion to terminate when `dt` drops below `√eps(T)`.

Pass `force_dtmin=true` to `solve` when using this callback.  Otherwise, the time-step size
may decrease too much *within* a single time step, so that the integrator itself will quit
before reaching this callback, leading to a less graceful exit.

If this terminator is triggered while `v` is less than 0.35, a warning will always be
issued; otherwise an `info` message will be issued only if the `quiet` flag is set to
`false`.
"""
function dtmin_terminator(PNType, T, quiet=false)
    sqrtϵ::T = √eps(T)  # Tricks for faster closures
    discrete_condition = let sqrtϵ = sqrtϵ
        (state, t, integrator) -> abs(integrator.dt) < sqrtϵ
    end
    function discrete_terminator!(integrator)
        v = integrator.u[symbol_index(PNType, :v)]
        message = (
            "Terminating evolution because the time-step size has become very small:\n" *
            "|dt=$(integrator.dt)| < √ϵ=$(sqrtϵ)\n" *
            "This is only unexpected for 𝑣 ≲ 0.35; the current value is 𝑣=$v."
        )
        if v < 7//20
            @warn message
        elseif !quiet
            @info message
        end
        return terminate!(integrator)
    end
    return DiscreteCallback(
        discrete_condition, discrete_terminator!; save_positions=(false, false)
    )
end

"""
    decreasing_v_terminator([quiet])

Construct termination criterion to stop integration when `v` is decreasing.

Note that some systems may truly have decreasing `v` as physical solutions — including
eccentric systems and possibly precessing systems.  You may prefer to implement another
solution, like detecting when `v` decreases below some threshold, or detecting when `v` is
decreasing too quickly.  See this function's source code for a simple

If this terminator is triggered while `v` is less than 0.35, a warning will always be
issued; otherwise an `info` message will be issued only if the `quiet` flag is set to
`false`.
"""
function decreasing_v_terminator(PNType, quiet=false)
    function discrete_condition(state, t, integrator)
        # How come here I can't just write get_du(integrator)[:v] ?
        return get_du(integrator)[symbol_index(PNType,:v)] < 0  # This translates to v̇<0
    end
    function discrete_terminator!(integrator)
        v = integrator.u[symbol_index(PNType,:v)]
        ∂ₜv = get_du(integrator)[symbol_index(PNType,:v)]
        message = (
            "Terminating forwards evolution because 𝑣 is decreasing:\n" *
            "This is only unusual if 𝑣 ≲ 0.35; the current value is 𝑣=$v\n" *
            "∂ₜ𝑣=$∂ₜv."
        )
        if v < 7//20
            @warn message
        elseif !quiet
            @info message
        end
        return terminate!(integrator)
    end
    return DiscreteCallback(
        discrete_condition, discrete_terminator!; save_positions=(false, false)
    )
end

"""
    nonfinite_terminator()

Construct termination criterion to stop integration when any NaN or Inf is found in the data
after an integration step.

If this terminator is triggered, a warning will always be issued.
"""
function nonfinite_terminator()
    function discrete_condition(state, t, integrator)
        return !(all(isfinite, state) && isfinite(t) && isfinite(integrator.dt))
    end
    function discrete_terminator!(integrator)
        @warn "Terminating forwards evolution because a non-finite number was found"
        return terminate!(integrator)
    end
    return DiscreteCallback(
        discrete_condition, discrete_terminator!; save_positions=(false, false)
    )
end
