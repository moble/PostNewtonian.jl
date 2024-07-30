# Computing the waveform

Once you have [computed the orbital evolution](@ref orbital_evolution), and
[selected the time steps](@ref sec_optional_choose_time_steps) on which you want
to evaluate the waveform, you can do so with one of the functions documented
below.

## Precise definition of waveforms

The waveform returned by these functions is essentially the complex strain ``h
\colonequals h_+ - i\,h_\times``.  However, this quantity decays as ``1/r``,
where ``r`` is the radius at which the strain is measured.  Therefore, as is
conventional, the returned quantity is actually

```math
H \colonequals \lim_{r\to\infty} h\, \frac{r}{M}\, \frac{c^2}{G},
```

where ``M`` is the total mass of the two objects in the binary.  Note that both
``h`` and ``H`` are dimensionless, but only ``H`` is scale invariant.

Binaries are usually modeled in relativity as being isolated systems in
otherwise empty asymptotically flat spacetimes — or even more specifically, in
asymptotically Minkowski spacetimes.  In this case, ``r`` is just the asymptotic
areal radius.  The equivalent expression in an FLRW spacetime requires a simple
reinterpretation of ``M`` as the redshifted mass ``M_z \colonequals M(1+z)`` and
``r`` as the luminosity distance ``d_L.``[^1]  Thus, to obtain the strain as it
would be measured (in the asymptotic approximation) by an actual observer in an
asymptotically flat universe or in our universe, we just need to invert the
previous equation:

```math
h \approx H\, \frac{M}{r}\, \frac{G}{c^2}
\qquad \mathrm{or} \qquad
h \approx H\, \frac{M_z}{d_L}\, \frac{G}{c^2}.
```

[^1]: Note that we use "luminosity distance" as it is [understood in
      contemporary cosmology](https://arxiv.org/abs/astro-ph/9905116), which is
      unfortunately different from its meaning in older cosmology literature,
      and even some of the gravitational-wave astronomy literature from before
      2015 or so.

Furthermore, ``H`` is a function of coordinates ``(t, \theta, \phi)``.  The
dependence on ``t`` will be discretely sampled at the times ``t_i`` that are
present in the `inspiral` arguments to the functions below.  To handle the
angular dependence, we provide the waveform decomposed as mode weights in a
spin-weighted spherical-harmonic decomposition, so that the actual quantity
returned will be

```math
H_{\ell,m}(t_i)
\colonequals
\int H(t, \theta, \phi)\, {}_{-2}\bar{Y}_{\ell,m}(\theta, \phi)\,
\sin\theta\, d\theta\, d\phi.
```

The output array is a two-dimensional complex array.  The first dimension varies
over ``(\ell,m)`` values, with ``m`` varying most rapidly — as in

```julia
[(ℓ,m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]
```

The second dimension of the array indexes the time.

See [Example 2 on the "Units" page](@ref Units-example-2) for a complete example
of converting this package's output to standard units in our universe.

## Evaluation of waveforms

```@docs
coorbital_waveform
inertial_waveform
```

## In-place evaluation of waveforms

This is likely to be an uncommon scenario, but if you happen to need to evaluate
the waveform repeatedly on the same set (or at least same number) of time steps
and same range of ``\ell`` values, then you can pre-allocate the storage needed
in this process, and re-use that storage.  To do so, pre-allocate using either
[`coorbital_waveform_computation_storage`](@ref) or
[`inertial_waveform_computation_storage`](@ref), and then pass the result as the
first argument to either [`coorbital_waveform!`](@ref) or
[`inertial_waveform!`](@ref).

!!! note
    The `h` array returned by either `coorbital_waveform!` or
    `inertial_waveform!` is part of the pre-allocated storage.  You need to use
    its values *before* you call either of those functions again, or those
    values will just be overwritten.

```@docs
coorbital_waveform_computation_storage
inertial_waveform_computation_storage
coorbital_waveform!
inertial_waveform!
```
