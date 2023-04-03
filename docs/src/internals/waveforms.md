## Computing the waveform from inspiral data

Once you have used [`orbital_evolution`](@ref) to compute the inspiral, and
selected the time steps on which you want the waveform, you can compute the
waveform with one of the following functions.

!!! danger
    Remember that, by default, [`orbital_evolution`](@ref) will return the
    solution at whatever time steps the ODE integrator chose, which are almost
    certainly too far apart to satisfy the Nyquist limit relevant for waveforms.
    You can use the `saves_per_orbit` argument to that function with a value of
    at least `2ℓₘₐₓ` (where `ℓₘₐₓ` is the argument to one of the functions
    below) to barely satisfy the Nyquist limit.  I usually use `4ℓₘₐₓ` for a
    little extra safety factor.

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
