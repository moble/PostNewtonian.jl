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
then be given in units of that mass.

The complete set of (required or optional) dimensionful arguments to
`orbital_evolution` is
```julia
M₁, M₂, Ωᵢ, λ₁, λ₂, Ω₁, Ωₑ, abstol, saveat
```
If we scale the values by some `σ`, then to maintain the same physical meaning
we should transform all the arguments as
```
M₁ ↦ M₁ * σ
M₂ ↦ M₂ * σ
Ωᵢ ↦ Ωᵢ / σ
λ₁ ↦ λ₁ * σ^5
λ₂ ↦ λ₂ * σ^5
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
All other evolved variables are dimensionless, and so are unaffected by the scaling.


## Example 1: Scale invariance

It's important to check that the claims above are actually true.  Here, we put
them to the test with a very large scale factor.

```@example units1
using PostNewtonian, Plots  # hide
plotlyjs()  # hide
default(size=(800,480), linewidth=3, leg=:top, legendfontsize=11)  # hide
σ = 10^10  # Scale factor

M₁ = 0.4
M₂ = 0.6
χ⃗₁ = [0.0, 0.5, 0.8]
χ⃗₂ = [0.8, 0.0, 0.5]
Ωᵢ = 0.01
Ω₁ = 3Ωᵢ/4
Ωₑ = 0.9
dt = 5

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
<iframe src="units1.html" style="height:500px;width:100%;"></iframe>
```

The fact that these curves are on top of each other, and the final expression
returned `true` shows that this system is indeed scale invariant.

## [Example 2: Astrophysical units](@id Units-example-2)

Suppose we want to construct a system with masses ``38.9\,M_\odot`` and
``32.7\,M_\odot``, as it would be observed from a distance of
``440\,\mathrm{Mpc}`` with time measured in seconds and dimensionless strain,
for an initial frequency in the ``(\ell,m)=(2,2)`` mode of ``20\,\mathrm{Hz}``.
(Note that these parameters are consistent with the estimated parameters of the
[GW150914](https://en.wikipedia.org/wiki/First_observation_of_gravitational_waves)
detection.)

Because the masses and frequencies must be in inverse units, we arbitrarily
choose to measure frequencies in their natural unit of Hertz, and therefore
masses are measured in seconds — multiplying the masses by ``G/c^3`` for
geometric units.  The distance to the source will also be converted to seconds
by dividing by ``c``, so that we can simply multiply the waveform by `(M₁+M₂)/r`
to get the observed strain.

```@example units2
using Quaternionic
using SphericalFunctions
using PostNewtonian
using Plots
plotlyjs()  # hide
default(size=(800,480), linewidth=3, leg=:top, legendfontsize=11)  # hide

# Useful astronomical constants
c = float(299_792_458) # m/s
GMₛᵤₙ = 1.32712440041e20 # m^3/s^2
au = float(149_597_870_700) # m
pc = 1au / (π / 648_000) # m
Mpc = 1_000_000pc # m

# Parameters of this system
M₁ = 38.9GMₛᵤₙ / c^3 # s
M₂ = 32.7GMₛᵤₙ / c^3 # s
r = 440Mpc / c # s
fᵢ = 20 # Hz

# Approximate maximum a posteriori spins
χ⃗₁ = [0.0, 0.0, 0.3224]
χ⃗₂ = [0.2663, 0.2134, -0.5761]

Ωᵢ = 2π * fᵢ / 2 # Hz — Note the factor of 2 because fᵢ represents the (2,2) mode
dt = 1/2048 # s — This is the sampling rate

inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ, saveat=dt)
hₗₘ = inertial_waveform(inspiral) * (M₁+M₂) / r

# Evaluate waveform at a point; see Scri.jl for simpler methods
R = Quaternionic.from_spherical_coordinates(2.856, 0.0)
₋₂Yₗₘ = SphericalFunctions.ₛ𝐘(-2, 8, Float64, [R])[1, :]
h = (₋₂Yₗₘ' * hₗₘ)[1,:]

plot(inspiral.t, real.(h), label="ℎ₊")
plot!(inspiral.t, -imag.(h), label="ℎₓ")
plot!(xlabel="Time (seconds)", ylabel="Strain (dimensionless)", ylim=(-1.5e-21,1.5e-21))
savefig("units2.html"); nothing  # hide
```
```@raw html
<iframe src="units2.html" style="height:500px;width:100%;"></iframe>
```
