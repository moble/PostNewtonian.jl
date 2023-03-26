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
    function conditions(out,state,t,integrator)
        out[1] = state[1]  # Terminate if M₁≤0
        out[2] = state[2]  # Terminate if M₂≤0
        out[3] = 1 - abs2vec(QuatVec{typeof(vₑ)}(state[3:5]...))  # Terminate if χ₁>1
        out[4] = 1 - abs2vec(QuatVec{typeof(vₑ)}(state[6:8]...))  # Terminate if χ₂>1
        out[5] = vₑ - state[13]  # Terminate at v = vₑ
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
    function terminators_backwards(out,state,t,integrator)
        out[1] = state[1]  # Terminate if M₁≤0
        out[2] = state[2]  # Terminate if M₂≤0
        out[3] = 1 - abs2vec(QuatVec{typeof(v₁)}(state[3:5]...))  # Terminate if χ₁>1
        out[4] = 1 - abs2vec(QuatVec{typeof(v₁)}(state[6:8]...))  # Terminate if χ₂>1
        out[5] = v₁ - state[13]  # Terminate at v = v₁
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
    function discrete_condition(state,t,integrator)
        abs(integrator.dt) < ϵ
    end
    function discrete_terminator!(integrator)
        @warn (
            "Terminating forwards evolution because time-step size is too small:\n"
            * "|dt=$(integrator.dt)| < √ϵ=$(ϵ)\n"
            * "This is probably fine if `v` ≳ 1/2."
        )
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
    function discrete_condition(state,t,integrator)
        # any(isnan, state) || isnan(t) || isnan(integrator.dt)
        !(all(isfinite, state) && isfinite(t) && isfinite(integrator.dt))
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
