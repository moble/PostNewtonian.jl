```@meta
CurrentModule = PostNewtonian
```

# PostNewtonian

This package computes orbital dynamics of and waveforms from binary black-hole
systems, in the [post-Newtonian
approximation](https://en.wikipedia.org/wiki/Post-Newtonian_expansion).
Currently, there is no support for eccentric systems, but general precessing
quasicircular systems are supported.

## Installation

This package requires Julia version 1.9 or greater.  If you haven't installed
Julia yet, you probably want to use
[`juliaup`](https://github.com/JuliaLang/juliaup#readme) to do so.  You'll
probably also want to use a Julia "project" specifically for this package.  An
easy way to do this is to create a directory, `cd` into that directory, and then run julia as
```
julia --project=.
```
Then, installation of this package involves the usual commands:
```julia
using Pkg
Pkg.add("PostNewtonian")
```

## Quick start

We can integrate the orbital dynamics of a black-hole binary using the
[`orbital_evolution`](@ref) function.  (Note that you don't have to use cool
Unicode names for your variables if you don't want to; `chi1` works just as
well as `χ⃗₁`, for example.)
```@example 1
using PostNewtonian

# Initial values of the masses and spins
M₁ = 0.6
M₂ = 0.4
χ⃗₁ = [0.1, 0.2, 0.3]
χ⃗₂ = [-0.3, -0.1, 0.7]

# Orbital angular frequency at the initial instant
Ωᵢ = 0.005

inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ);
nothing  # hide
```
There are numerous keyword arguments to `orbital_evolution`.  You may wish to
pass the `quiet=true` option, for example, to silence informational messages
like the one above about successful evolutions.  Other keyword arguments control
things like the range of frequencies over which to integrate (including possibly
both backwards *and* forwards from `Ωᵢ`), accuracy of the ODE integration, time
steps at which to save the results, the PN order at which to compute the
evolution equations, and so on.  See that function's [documentation](@ref
orbital_evolution) for details.

The returned `inspiral` object is an
[`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) object
with many features, like high-order interpolation (using dense output from the
ODE integrator).  The time steps at which the solution was saved are available
as `inspiral.t`, and the evolved variables are available as `inspiral.u`, or by
their names as in `inspiral[:v]` or `inspiral[:Φ]`.  For example, we can plot
the components of the spin of object 1 like this:
```@example 1
using Plots
plotlyjs()  # Requires also adding `PlotlyJS` to your project

plot(
    inspiral, idxs=[(0,:χ⃗₁ˣ), (0,:χ⃗₁ʸ), (0,:χ⃗₁ᶻ)], leg=:top, linewidth=3,
    xlabel="Time (M)", ylabel="Dimensionless spin components"
)
savefig("inspiral_spins.html"); nothing  # hide
```
```@raw html
<iframe src="inspiral_spins.html" style="height:500px;width:100%;"></iframe>
```

Usually, we will also want the actual waveform from this system.  Using the
`inspiral` object from above, we just call
```julia
# @example 1
h = mode_weights(inspiral)
plot(
    h, idxs=[(0,:χ⃗₁ˣ), (0,:χ⃗₁ʸ), (0,:χ⃗₁ᶻ)], leg=:top, linewidth=3,
    xlabel="Time (M)", ylabel="Mode weights"
)
savefig("mode_weights.html"); nothing  # hide
```
```
@raw html
<iframe src="mode_weights.html" style="height:500px;width:100%;"></iframe>
```
