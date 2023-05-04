# High-level interface

The typical workflow for users of this package will involve two or three steps:

1. Integrate the orbital dynamics of a black-hole binary using the
   [`orbital_evolution`](@ref) function.

2. (Optional) If not specified in step 1, choose the time steps on which you
   want the waveform.

3. Compute the waveform as a function of the orbital dynamics using the
   [`inertial_waveform`](@ref) (or [`coorbital_waveform`](@ref)) function.

Here, we'll work through an example, and provide some more details.

## 1. Integrate orbital dynamics

First, we have to specify the initial masses, spins, and orbital angular
frequency.  Let's arbitrarily choose something close to a [hangup-kick](@ref
hangup_kick) configuration.
```@example 1
using PostNewtonian

# Initial values of the masses, spins, and orbital angular frequency
M₁ = 0.6
M₂ = 0.4
χ⃗₁ = [0.7, 0.1, 0.7]
χ⃗₂ = [-0.7, 0.1, 0.7]
Ωᵢ = 0.01

inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ)
nothing;  # hide
```
There are also numerous optional keyword arguments to `orbital_evolution`,
controlling things like the range of frequencies over which to integrate
(including possibly both forwards *and* backwards from the initial values),
accuracy of the ODE integration, the PN order at which to compute the evolution
equations, and so on.  See that function's [documentation](@ref
orbital_evolution) for details.

The returned object named `inspiral` here is an
[`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) with
many features, like high-order interpolation (using dense output from the ODE
integrator).  The time steps at which the solution was saved are available as
`inspiral.t`, and the evolved variables are available as `inspiral.u`, or by
their names as in `inspiral[:v]` or `inspiral[:Φ]`.  For example, we can plot
the components of the spin of object 1 like this:
```@example 1
using Plots  # Requires also installing `Plots` in your project
plotlyjs()  # hide
default(size=(800,480), linewidth=2, leg=:top)  # hide

plot(inspiral.t, inspiral[:χ⃗₁ˣ], label=raw"$\vec{\chi}_1^x$")
plot!(inspiral.t, inspiral[:χ⃗₁ʸ], label=raw"$\vec{\chi}_1^y$")
plot!(inspiral.t, inspiral[:χ⃗₁ᶻ], label=raw"$\vec{\chi}_1^z$")
plot!(xlabel="Time (𝑀)", ylabel="Dimensionless spin components")
savefig("inspiral_spins.html"); nothing  # hide
```
```@raw html
<iframe src="../inspiral_spins.html" style="height:500px;width:100%;"></iframe>
```
As expected, we see *significant* precession of the spin on long time scales, as
well as smaller nutations on orbital time scales visible mostly at later times.

The evolved variables, in order, are
```@example 1
join(stdout, PostNewtonian.pnsystem_symbols, ", ")  # hide
```
They can be accessed by their symbols, like the spins above, or by their number
in this list.  To access the `i`th variable at time step `j`, use `sol[i, j]`.
You can also use colons: `sol[i, :]` is a vector of the `i`th variable at all
times, and `sol[:, j]` is a vector of all the data at time step `j`.  For
example, `inspiral[:χ⃗₁ˣ]` could also be written as `inspiral[3, :]`.

By default, the output of `orbital_evolution` is just the time steps to which
the adaptive ODE integrator happened to step.  If you know that you want the
solution on a set of uniform time steps separated by `dt` — such as when you
need to FFT the waveform — you can pass the option `saveat=dt`.  Or, if you
somehow know the specific set of times `t` at which you want the solution, you
can pass `saveat=t`.  Finally, if you want the solution — say — 32 times per
orbit, you can pass the option `saves_per_orbit=32`, which calls
[`uniform_in_phase`](@ref) as needed.

## 2. (Optional) Choose time steps

If you did not pass the `saveat` or `saves_per_orbit` arguments to
`orbital_evolution` (as described in the previous paragraph), the output will
usually be on a fairly coarse set of times — possibly many times the orbital
period for non-precessing systems.  However, in this case the solution will come
with "dense output", which lets us quickly and accurately interpolate the
solution to a new set of times.  This is important if you want the waveform to
be sampled at the "local" [Nyquist
rate](https://en.wikipedia.org/wiki/Nyquist_rate); anything less and the
waveform will be distorted by
[aliasing](https://en.wikipedia.org/wiki/Aliasing).

For instance, suppose we want to plot the results every ``0.5M`` for the last
``5,000M`` of the data.  We could define the new set of times as
```@example 1
t′ = inspiral.t[end]-5_000 : 0.5 : inspiral.t[end]
nothing;  # hide
```
and then interpolate to this new set of times as
```@example 1
inspiral = inspiral(t′)
nothing  # hide
```
We couldn't have achieved quite the same effect with the `saveat` argument
mentioned above because in this case, we wanted to know the point at which we
were ``5,000M`` before the end of the inspiral.  In general, if you need the
results of `orbital_evolution` before you can decide what you want `t′` to be,
this is the approach you'll have to take.

## 3. Compute the waveform

Usually, we will also want the actual waveform from this system.  We can just
call [`inertial_waveform`](@ref) (or [`coorbital_waveform`](@ref) for the
waveform in a rotating frame in which the binary is stationary).
```@example 1
h = inertial_waveform(inspiral)
nothing;  # hide
```
Again, we can plot the result:
```@example 1
plot(inspiral.t, real.(h[1, :]), label=raw"$\Re\left\{h_{2,2}\right\}$")
plot!(inspiral.t, imag.(h[1, :]), label=raw"$\Im\left\{h_{2,2}\right\}$")
plot!(inspiral.t, abs.(h[1, :]), label=raw"$\left|h_{2,2}\right|$", linewidth=3)
plot!(xlabel="Time (𝑀)", ylabel="Mode weights", ylim=(-0.5,0.5))
savefig("waveform.html"); nothing  # hide
```
```@raw html
<iframe src="../waveform.html" style="height:500px;width:100%;"></iframe>
```
We see oscillations in the amplitude of the ``h_{2,2}`` mode on the orbital
timescale, which is to be expected in a hangup-kick scenario as the system
alternates between beaming power preferentially along the ``+z`` and ``-z``
directions.
