```@meta
CurrentModule = PostNewtonian
```

# PostNewtonian

This package computes orbital dynamics of and waveforms from binary black-hole
systems, in the [post-Newtonian
approximation](https://en.wikipedia.org/wiki/Post-Newtonian_expansion).
Currently, general precessing quasispherical systems are supported, but support
for eccentric systems is still upcoming.

## Installation

If you intend to use this package via Python, see [this page](@ref
Using-this-package-from-Python) for installation instructions.

It is recommended to use this package with Julia version 1.9 or greater, because
of that version's improved pre-compilation caching.  If you find it very slow
the first time you use functions from this package in a new Julia session, that
is most likely because Julia has to compile a lot of code.  Version 1.9 does a
better job of caching that compiled code, which speeds up your first-time usage.

If you haven't installed Julia yet, you probably want to use
[`juliaup`](https://github.com/JuliaLang/juliaup#readme) to do so.  You'll
probably also want to use a Julia ["project
environment"](https://pkgdocs.julialang.org/v1/environments/) specifically for
using this package.  An easy way to do this is to create a directory, `cd` into
that directory, and then run julia as
```
julia --project=.
```
Then, installation of this package involves the usual commands:
```julia
using Pkg
Pkg.add("PostNewtonian")
```


## Quick start

An example with slightly more explanation is given under ["High-level
interface"](@ref High-level-interface), and of course the rest of this
documentation goes into far more detail.  Here we see a simple example to start
things off.

!!! tip
    You don't have to use [cool Unicode
    names](https://docs.julialang.org/en/v1/manual/unicode-input/) for
    your variables if you don't want to.  For example, `chi1` works just
    as well as `χ⃗₁`.  Similarly, many functions in this package have
    Unicode names or take optional Unicode keyword arguments.  But every
    such name or argument will also have an ASCII equivalent; see the
    documentation of those functions for the appropriate substitutions.

```@example 1
using PostNewtonian

# Initial values of the masses, spins, and orbital angular frequency
M₁ = 0.4
M₂ = 0.6
χ⃗₁ = [0.0, 0.5, 0.8]
χ⃗₂ = [0.8, 0.0, 0.5]
Ωᵢ = 0.01

# Integrate the orbital dynamics
inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ)

# Interpolate for nicer plotting
t′ = inspiral.t[end]-6_000 : 0.5 : inspiral.t[end]
inspiral = inspiral(t′)

# Compute the waveform in the inertial frame
h = inertial_waveform(inspiral)
nothing  # hide
```
We can plot the result like this:
```@example 1
using Plots  # Requires also installing `Plots` in your project
plotlyjs()  # hide
default(size=(800,480), linewidth=2, leg=:top, legendfontsize=11)  # hide
default(extra_plot_kwargs = KW(:include_mathjax => "cdn"))  # hide

plot(inspiral.t, real.(h[1, :]), label=raw"$\Re\left\{h_{2,-2}\right\}$")
plot!(inspiral.t, imag.(h[1, :]), label=raw"$\Im\left\{h_{2,-2}\right\}$")
plot!(inspiral.t, abs.(h[1, :]), label=raw"$\left|h_{2,-2}\right|$", linewidth=3)
plot!(inspiral.t, abs.(h[5, :]), label=raw"$\left|h_{2,2}\right|$")
plot!(xlabel=raw"$\text{Time }(M)$", ylabel="Mode weights", ylim=(-0.45,0.45))

savefig("quickstart_waveform.html"); nothing  # hide
```
```@raw html
<iframe src="quickstart_waveform.html" style="height:500px;width:100%;"></iframe>
```
We see various features to be expected of a precessing system like this,
including slow modulations of the modes on the precession timescale, as well as
faster oscillations in the amplitudes and asymmetry between the ``m=\pm 2``
modes on the orbital timescale.
