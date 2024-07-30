# A note about units

!!! note "TL;DR:"
    This package uses geometric units with ``G=c=1``.  Units are essentially
    determined by the values you input.  Masses and frequencies in particular
    must have inverse units.  For example, if frequency is measured in Hertz,
    mass should be measured in seconds.  The common relativist's convention of
    measuring in units of "total mass" (whatever that means to you) is also
    perfectly fine.

    Output time is in units of the total mass; to get the time in seconds,
    multiply by the total mass as measured in seconds.  Waveforms are output
    as the asymptotic quantity ``rh/M``; to get an observable strain,
    multiply by the total mass and divide by the distance to the source.

    You may want to skip to [Example 2](@ref Units-example-2) for a complete
    example using astrophysical units.

The units of measurement used in this package are implicitly established by the
input.  Specifically, the arguments `M₁`, `M₂`, and `Ωᵢ` in the call

```julia
inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ)
```

are inherently dimensionful, and the units are effectively determined by the
values entered.  For example, if we enter `0.4` as `M₁`, we have established the
units as — quite simply — those units in which `M₁` has the value `0.4`.  The
only real requirement on the units of the input arguments is that they must be
consistent.  In particular, the quantity `(M₁+M₂)*Ωᵢ` should be dimensionless.

It must also be noted that this package uses geometric units where ``G=c=1``; in
fact, ``G`` and ``c`` never appear in the code at all.  This means that you must
input values in geometric units, and interpret the output as being in geometric
units.  This does not completely constrain the units, however.  There is still a
freedom to scale the units however you wish.  For example, it is most common to
describe all quantities in units of some mass ``M`` — such as the sum of the
input masses `M₁+M₂` or the Solar mass ``M_\odot``.  The input arguments would
then *all* be given in units of that mass, or more likely some product of that
mass with ``G`` and ``c``.

The complete set of (required or optional) dimensionful arguments to
`orbital_evolution` is

```julia
M₁, M₂, Ωᵢ, Λ₁, Λ₂, Ω₁, Ωₑ, abstol, saveat
```

If we scale the values by some `σ`, then to maintain the same physical meaning
we should transform all the arguments as

```julia
M₁ ↦ M₁ * σ
M₂ ↦ M₂ * σ
Ωᵢ ↦ Ωᵢ / σ
Λ₁ ↦ Λ₁
Λ₂ ↦ Λ₂
Ω₁ ↦ Ω₁ / σ
Ωₑ ↦ Ωₑ / σ
abstol ↦ [abstol[1:2] * σ; abstol[3:end]]
saveat ↦ saveat * σ
```

Note that the scaling happens automatically for default values of the
parameters; you would only need to rescale them if you were actually setting
them.

When the scaled arguments are provided to `orbital_evolution`, the resulting
`inspiral` is affected as

```julia
inspiral.t ↦ inspiral.t * σ
inspiral[:M₁] ↦ inspiral[:M₁] * σ
inspiral[:M₂] ↦ inspiral[:M₂] * σ
```

All other evolved variables are dimensionless, and so are unaffected by the
scaling.  In particular, if the initial conditions are entered with values of
`M₁` and `M₂` in their desired final units, no change needs to be made.  On the
other hand, if the values are entered so that `M₁+M₂=1`, we can rescale any BBH
system to any desired total mass `Mₜₒₜ` using

```julia
inspiral.t ↦ inspiral.t * Mₜₒₜ * G / c^3
```

Furthermore, if `M₁` and `M₂` represent the *source-frame* mass, and we need to
apply a cosmological redshift `z`, we simply apply

```julia
inspiral.t ↦ inspiral.t * (1+z)
```

The waveform ``H`` returned by either [`coorbital_waveform`](@ref) or
[`inertial_waveform`](@ref) is the rescaled asymptotic limit of the strain.  The
observable strain is

```math
h \approx H\, \frac{M}{r}\, \frac{G}{c^2}
\qquad \mathrm{or} \qquad
h \approx H\, \frac{M_z}{d_\mathrm{L}}\, \frac{G}{c^2}.
```

(The equality is only approximate because we assume that terms in ``1/r^2`` can
be ignored, which is almost certainly true of any detection in the foreseeable
future.)  Here, ``r`` is the radius of the observer in asymptotically Minkowski
coordinates centered on the source and ``M`` is the total mass of the binary;
alternatively, ``M_z=M(1+z)`` is the redshifted mass and ``d_\mathrm{L}`` is the
luminosity distance between the source and observer.  For more complete
description, see [here](@ref Computing-the-waveform).

## Example 1: Scale dependence

It's important to check that the claims above are actually true.  Here, we put
them to the test with a very large scale factor.

```@example units1
using PostNewtonian, Plots  # hide
plotlyjs()  # hide
default(size=(800,480), linewidth=3, leg=:top, legendfontsize=11)  # hide
σ = 10.0^10  # Scale factor

M₁ = 0.4
M₂ = 0.6
χ⃗₁ = [0.0, 0.5, 0.8]
χ⃗₂ = [0.8, 0.0, 0.5]
Ωᵢ = 0.01
Ω₁ = 3Ωᵢ/4
Ωₑ = 0.9
dt = 5.0

inspiral1 = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, Ω₁=Ω₁, Ωₑ=Ωₑ, saveat=dt)
inspiral2 = orbital_evolution(M₁*σ, M₂*σ, χ⃗₁, χ⃗₂, Ωᵢ/σ, Ω₁=Ω₁/σ, Ωₑ=Ωₑ/σ, saveat=dt*σ)

plot(inspiral1.t, inspiral1[:v], label="Original")
plot!(inspiral2.t/σ, inspiral2[:v], label="Scaled", linewidth=4, ls=:dot)
plot!(xlabel="Time", ylabel="PN parameter 𝑣")
savefig("units1.html"); nothing  # hide

# Check that the evolved `v` parameters are nearly equal (when interpolating
# the second to *exactly* the same times as the first)
inspiral1[:v] ≈ inspiral2(inspiral1.t*σ, idxs=:v)
```

```@raw html
<!-- NOTE: ../ in src works on github, but not locally -->
<iframe src="../units1.html" style="height:500px;width:100%;"></iframe>
```

The fact that these curves are on top of each other, and the final expression
returned `true` shows that this system does indeed depend on scale as suggested
above.

## [Example 2: Astrophysical units](@id Units-example-2)

Suppose we want to construct a system consistent with
[GW150914](https://arxiv.org/abs/1606.01210), having source-frame masses
``35.6\,M_\odot`` and ``30.0\,M_\odot``, as it would be observed from a distance
of ``440\,\mathrm{Mpc}`` with an initial frequency in the dominant
``(\ell,m)=(2,\pm 2)`` modes of ``f_i \approx 20\,\mathrm{Hz}``.

Because the masses and frequencies must be in inverse units, we arbitrarily
choose to measure frequencies in their natural unit of Hertz, and therefore
masses are measured in seconds — multiplying the masses by ``G/c^3`` for
geometric units.  The distance to the source ``d_\mathrm{L}`` will also be
converted to seconds by dividing by ``c``, so that we can simply multiply the
waveform by ``M_z/d_\mathrm{L}`` to get the observed strain.

The frequency ``f_i`` is the observed initial frequency of the ``(2,2)`` mode.
We will need the corresponding angular orbital frequency in the frame.  Recall
that, by definition of redshift ``z``, we have

```math
\frac{f_{\mathrm{source}}} {f_{\mathrm{observer}}} = 1+z
\qquad
\text{and}
\qquad
\frac{\Delta t_{\mathrm{observer}}}{\Delta t_{\mathrm{source}}} = 1+z.
```

Thus, the initial *source-frame orbital angular frequency* ``\Omega_i`` is
related to the *observer-frame $(2,2)$ temporal frequency* ``f_i`` by[^1]

```math
%\Omega_i = 2\pi f_{\mathrm{source}}^{(2,2)} / 2 \qquad f_i = f_{\mathrm{observer}}^{(2,2)}
%\qquad
%\frac{f_{\mathrm{source}}^{(2,2)}} {f_{\mathrm{observer}}^{(2,2)}} =
%\frac{\Omega_i/\pi} {f_i} = 1+z
%\qquad
\Omega_i = \pi\, f_i\, (1+z).
```

[^1]: Here, we have approximated the frequency of the ``(\ell,m)=(2,\pm 2)``
    modes as being twice the orbital frequency.  In the post-Newtonian
    approximation, waveform mode ``H_{\ell,m}`` is frequently written as
    ``A_{\ell,m}\, \exp[-i\, m\, \Phi]``, where ``\Phi`` is the orbital phase
    that varies in time with ``d\Phi/dt = \Omega``.  Just looking at the
    exponential, we might expect the frequency of ``H_{\ell,m}`` to be ``m\,
    \Omega``.  However, the "amplitude" term ``A_{\ell,m}`` is not an amplitude
    in the usual sense; it is actually complex, with a slowly varying phase.
    This means that the frequency of ``H_{\ell,m}`` is not precisely ``m\,
    \Omega`` — though this is a reasonable approximation for systems with
    reasonably low frequency.

```@example units2
using Quaternionic
using SphericalFunctions
using PostNewtonian
using Plots
plotlyjs()  # hide
default(size=(800,480), linewidth=3, leg=:top, legendfontsize=11)  # hide

# Useful astronomical constants
c = float(299_792_458)   # meters/second
GMₛᵤₙ = 1.32712440041e20  # meters³/second²
au = float(149_597_870_700)  # meters
pc = 1au / (π / 648_000)  # meters
Mpc = 1_000_000pc  # meters

# Approximate maximum a posteriori spins (via SXS:BBH:0308)
χ⃗₁ = [0.0, 0.0, 0.3224]
χ⃗₂ = [0.2663, 0.2134, -0.5761]

# Parameters of this system (arXiv:1606.01210)
M₁ = 35.6GMₛᵤₙ / c^3  # seconds
M₂ = 30.0GMₛᵤₙ / c^3  # seconds
dL = 440Mpc / c  # seconds
z = 0.094
fᵢ = 20  # Hertz
Δtₒ = 1/2048  # seconds

# Conversions
Mz = (M₁+M₂) * (1+z)  # seconds
Ωᵢ = π * fᵢ * (1+z)  # Hertz
Δtₛ = Δtₒ / (1+z)  # seconds

# Evaluate waveform and convert to observer's frame and units
inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, saveat=Δtₛ)
tₒ = (1+z) * inspiral.t
Hₗₘ = inertial_waveform(inspiral)
hₗₘ = Hₗₘ * Mz / dL

# Evaluate waveform at a point; see Scri.jl for simpler methods
R = Quaternionic.from_spherical_coordinates(2.856, 0.0)
₋₂Yₗₘ = SphericalFunctions.ₛ𝐘(-2, 8, Float64, [R])[1, :]
h = (₋₂Yₗₘ' * hₗₘ)[1,:]

# Plot the results
plot(tₒ, real.(h), label="ℎ₊")
plot!(tₒ, -imag.(h), label="ℎₓ")
plot!(xlabel="Time (seconds)", ylabel="Strain (dimensionless)", ylim=(-1.5e-21,1.5e-21))
plot!(title="Observer-frame waveform", legend=(0.205, 0.9))
savefig("units2.html"); nothing  # hide
```

```@raw html
<!-- NOTE: ../ in src works on github, but not locally -->
<iframe src="../units2.html" style="height:500px;width:100%;"></iframe>
```
