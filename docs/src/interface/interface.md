# High-level functions

Most users will just want to evaluate the PN orbital evolution and/or the
waveform.  Though these were both touched on briefly in the "[Quick
start](@ref)" section and are extensively documented in their corresponding
docstrings, we highlight some of the more useful features here.

## Orbital evolution

A post-Newtonian inspiral is computed by the [`orbital_evolution`](@ref)
function.

### [Input arguments](@id Orbital_evolution_Input_arguments)

#### Time steps

The default output contains just the times to which the ODE integrator stepped,
which can be *very* coarsely sampled.  In particular, if you evaluate the
inertial-frame waveform on these time steps, it will almost surely be sampled
far below the Nyquist limit, which will generally render the waveform useless.

There are three useful approaches to deal with this:

  1. Do nothing on input, but interpolate the *output* inspiral object to some
     desired time steps using its built-in interpolation capabilities.  This
     uses the "dense output" results of the ODE integrator, which produces very
     efficient and accurate results.  This is a very good option if the
     following two are not applicable, especially if you need to know the PN
     evolution before you can decide on the time steps you want.  See the
     ["Output"](@ref Orbital_evolution_Output) section below for details.
  2. If you want to ensure that each orbit is sampled a certain number of times,
     you can pass the `saves_per_orbit` option to `orbital_evolution`.
  3. If you want the results with some fixed uniform time step `dt`, you can
     just pass the option `saveat=dt` to `orbital_evolution`.
  4. If you happen to know specific times `t` at which you want the result, you
     can pass the option `saveat=t` to `orbital_evolution`.  This may be
     relevant, for example, when comparing the PN waveform to an existing
     waveform.  Note that the inspiral *may* end earlier than your last `t`
     value, in which case the output will just be shorter than requested.

#### When to stop

Once the equations have been entered correctly, the hardest part about PN
orbital evolution is deciding when to stop.  We know that PN approximations fail
somewhere below ``v=1``.  In practice, especially at higher orders, the
equations *frequently* for ``v \approx 0.5``, and can fail as early as ``v
\approx 0.4``.  The breakdown point typically depends on many factors, such as
the masses and spins of the binary, the PN order at which expressions are
truncated, and the approximant used to integrate the orbit.

But wherever this breakdown lies, as the ODE integrator approaches the
breakdown, it works harder to maintain accurate evolution of the equations it
has been given — even though *the equations themselves* are no longer
trustworthy.  In fact, a majority of the ODE integrator's time may be spent
trying to evolve the last fraction of an orbit as the integrator takes smaller
and smaller time steps to maintain accuracy in a regime where we shouldn't
believe PN itself.  Loosening the `abstol` and/or `reltol` requirements would
make the rest of the inspiral less accurate, but still wouldn't improve this
situation very much.

Ideally, you should pick some frequency below the default value (which
corresponds to ``v=1``) to let the ODE integration stop before it gets too close
to the breakdown point.  You can pass this lower value as the `Ωₑ` argument to
`orbital_evolution`.  For example, a value corresponding to ``v=0.3`` would be
unlikely to run into problems — except perhaps in very obscure corners of
parameter space.



### [Output](@id Orbital_evolution_Output)

`inspiral`
* inspiral.t
* inspiral[:, j]
* inspiral[i, j]
* inspiral[i, :]
* inspiral(t′)
* inspiral(t′, Val{1})  # Up to Val{4}
* inspiral(t′, Val{1}, idxs=[:M₁, :M₂])
* inspiral[:v]
* inspiral[[:v, :Φ], :]
  Note that there is no translation for this in Python (because of Unicode and because of the Symbol syntax), but we can use the following from Python:
* PNVariable.(inspiral.u)

## Waveform

`coorbital_waveform` / `inertial_waveform`

### [Input arguments](@id Waveform_Input_arguments)

### [Output](@id Waveform_Output)

