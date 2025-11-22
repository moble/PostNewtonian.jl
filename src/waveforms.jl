# Helper functions to get the `pnsystem` from either an ODESolution or the result of
# interpolating such a thing
_pnsystem(inspiral::ODESolution) = inspiral.prob.p
_pnsystem(inspiral::DiffEqArray) = inspiral.p

"""
    coorbital_waveform_computation_storage(inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    coorbital_waveform_computation_storage(inspiral; [ell_min=2], [ell_max=8], [PNOrder])

Construct storage needed to compute waveforms in the co-orbital frame without allocations.

This returns the storage for the waveforms themselves and `PNSystem` used as temporary
storage.  The returned quantity can just be passed as the first argument to
[`coorbital_waveform!`](@ref) without being unpacked.

The meaning of the arguments is the same as in [`coorbital_waveform`](@ref).
"""
function coorbital_waveform_computation_storage(
    inspiral;
    ell_min=2,
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(_pnsystem(inspiral)),
)
    @assert 0 ≤ ℓₘᵢₙ ≤ 2
    @assert ℓₘᵢₙ ≤ ℓₘₐₓ
    p = _pnsystem(inspiral)
    PNSystemType = parameterless_type(p)
    pnsystem = PNSystemType(copy(inspiral.u[1]); Λ₁=Λ₁(p), Λ₂=Λ₂(p), PNOrder)
    n_modes = (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2
    h = Array{Complex{eltype(inspiral)}}(undef, n_modes, length(inspiral))
    return h, pnsystem
end

# Helper function for default ℓₘᵢₙ associated to a certain modes_function
_default_ℓₘᵢₙ(::typeof(Ψ_M!)) = 0
_default_ℓₘᵢₙ(::typeof(h!)) = 2
# Don't have any other details... user must provide

"""
    coorbital_waveform!(storage, inspiral; [modes_function=h!], [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    coorbital_waveform!(storage, inspiral; [modes_function=h!], [ell_min=2], [ell_max=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights from the `modes_function` in
the co-orbital frame for the given `inspiral` output by
[`orbital_evolution`](@ref), using pre-allocated storage.

The storage is assumed to be the object returned from
[`coorbital_waveform_computation_storage`](@ref).  Other arguments are the same as in
[`coorbital_waveform`](@ref).
"""
function coorbital_waveform!(
    storage,
    inspiral;
    modes_function=(h!),
    ell_min=_default_ℓₘᵢₙ(modes_function),
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(storage[2]),
)
    h, pnsystem = storage
    @assert length(pnsystem.state) == length(inspiral.u[1])
    @assert length(inspiral) == size(h, 2)
    @assert (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2 == size(h, 1)
    @inbounds @fastmath for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        modes_function(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
    end
    return h
end

"""
    coorbital_waveform(inspiral; [modes_function=h!], [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    coorbital_waveform(inspiral; [modes_function=h!], [ell_min=2], [ell_max=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights from the `modes_function` in
the co-orbital frame for the given `inspiral` output by
[`orbital_evolution`](@ref).

See also [`inertial_waveform`](@ref) for the waveform in the inertial frame.

For a detailed description of the format and physical interpretation of the returned
quantity, see the ["Waveforms" section of the
manual.](https://moble.github.io/PostNewtonian.jl/stable/internals/waveforms.html)

The `PNOrder` defaults to the one used to compute `inspiral`, but may be changed by passing
the keyword argument.

!!! tip
    If you need this waveform at a different set of times `t′` than is currently present in
    `inspiral.t`, you should use the built-in interpolation capabilities of `inspiral`
    *first*, as in `inspiral′ = inspiral(t′)`, rather than interpolating the results of this
    function.  Or, perhaps better yet, you could select the times when calling
    `orbital_evolution` by using the `saves_per_orbit` or `saveat` keyword argument to that
    function.  These approaches will be more accurate, faster, and require less memory than
    interpolating the result of this function.  If using `saves_per_orbit`, you probably
    want to set it to *at least* `2ℓₘₐₓ`, or preferably `4ℓₘₐₓ`.

The mode weights are given starting at `ℓₘᵢₙ` (which must satisfy `0 ≤ ℓₘᵢₙ ≤ 2`) and
extending through `ℓₘₐₓ`.  The waveform is returned as a 2-dimensional `Array`, in which the
first index varies over the mode index from `(ℓ, m) = (ℓₘᵢₙ, -ℓₘᵢₙ)` to `(ℓ, m) = (ℓₘₐₓ,
ℓₘₐₓ)`, with `m` varying most rapidly, and the second index varying over the time steps.

In this function, the waveform is returned in the co-orbital frame — which is somewhat like
the co-rotating frame.  In particular, the modes of non-precessing systems vary slowly, over
inspiral timescales; modes of precessing systems still vary on orbital timescales, though
even this variation could be factored out.
"""
function coorbital_waveform(
    inspiral;
    modes_function=(h!),
    ell_min=_default_ℓₘᵢₙ(modes_function),
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(_pnsystem(inspiral)),
)
    storage = coorbital_waveform_computation_storage(inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
    return coorbital_waveform!(storage, inspiral; modes_function, ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
end

"""
    inertial_waveform_computation_storage(inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    inertial_waveform_computation_storage(inspiral; [ell_min=2], [ell_max=8], [PNOrder])

Construct storage needed to compute waveforms in the inertial frame without allocations.

This returns the storage for the waveforms themselves, storage used for computing the Wigner
𝔇 matrices, and for "in-place" multiplication.  The returned quantity can just be passed as
the first argument to [`inertial_waveform!`](@ref) without being unpacked.

The meaning of the arguments is the same as in [`inertial_waveform`](@ref).
"""
function inertial_waveform_computation_storage(
    inspiral;
    ell_min=2,
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(_pnsystem(inspiral)),
)
    h, pnsystem = coorbital_waveform_computation_storage(inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
    (D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ) = Dprep(ℓₘₐₓ, eltype(inspiral))
    hᵢ = similar(h[:, 1])
    return h, pnsystem, D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ, hᵢ
end

"""
    inertial_waveform!(storage, inspiral; [modes_function=h!], [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    inertial_waveform!(storage, inspiral; [modes_function=h!], [ell_min=2], [ell_max=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights from the `modes_function` in
the inertial frame for the given `inspiral` output by
[`orbital_evolution`](@ref), using pre-allocated storage.

The storage is assumed to be the object returned from
[`inertial_waveform_computation_storage`](@ref).  Other arguments are the same as in
[`inertial_waveform`](@ref).
"""
function inertial_waveform!(
    storage,
    inspiral;
    modes_function=(h!),
    ell_min=_default_ℓₘᵢₙ(modes_function),
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(storage[2]),
)
    h, pnsystem, D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ, hᵢ = storage
    @assert length(pnsystem.state) == length(inspiral.u[1])
    @assert length(inspiral) == size(h, 2)
    @assert (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2 == size(h, 1)
    @inbounds @fastmath for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        modes_function(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
        D!(D, conj(R(pnsystem)), ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
        f′ = Yiterator(hᵢ, ℓₘₐₓ, ℓₘᵢₙ, 1)
        f = Yiterator(h[:, iₜ], ℓₘₐₓ, ℓₘᵢₙ, 1)
        𝔇 = Diterator(D, ℓₘₐₓ, ℓₘᵢₙ; warn=false)
        for (f′ˡ, fˡ, 𝔇ˡ) ∈ zip(f′, f, 𝔇)
            mul!(f′ˡ, 𝔇ˡ, fˡ)
        end
        h[:, iₜ] .= hᵢ
    end
    return h
end

"""
    inertial_waveform(inspiral; [modes_function=h!], [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])
    inertial_waveform(inspiral; [modes_function=h!], [ell_min=2], [ell_max=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights from the `modes_function` in
the inertial frame for the given `inspiral` output by
[`orbital_evolution`](@ref).

The inertial frame is the one in which inertial observers are found, so this waveform is
more like one that actual observers would detect.  This function transforms the waveform
from the co-orbital frame — which is the one in which PN expressions are provided.

See [`coorbital_waveform`](@ref) for details about the other arguments and returned
quantity.
"""
function inertial_waveform(
    inspiral;
    modes_function=(h!),
    ell_min=_default_ℓₘᵢₙ(modes_function),
    ell_max=8,
    ℓₘᵢₙ=ell_min,
    ℓₘₐₓ=ell_max,
    PNOrder=pn_order(_pnsystem(inspiral)),
)
    storage = inertial_waveform_computation_storage(inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
    return inertial_waveform!(storage, inspiral; modes_function, ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
end
