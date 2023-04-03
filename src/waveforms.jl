# Helper functions to get the `pnsystem` from either an ODESolution or the result of
# interpolating such a thing
_pnsystem(inspiral::SciMLBase.ODESolution) = inspiral.prob.p
_pnsystem(inspiral::SciMLBase.DiffEqArray) = inspiral.p


"""
    coorbital_waveform_computation_storage(inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Construct storage needed to compute waveforms in the co-orbital frame without allocations.

This returns the storage for the waveforms themselves and `PNSystem` used as temporary
storage.  The returned quantity can just be passed as the first argument to
[`coorbital_waveform!`](@ref) without being unpacked.

The meaning of the arguments is the same as in [`coorbital_waveform`](@ref).
"""
function coorbital_waveform_computation_storage(
    inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(_pnsystem(inspiral))
)
    @assert 0 ≤ ℓₘᵢₙ ≤ 2
    @assert ℓₘᵢₙ ≤ ℓₘₐₓ
    p = _pnsystem(inspiral)
    PNSystemType = SciMLBase.parameterless_type(p)
    pnsystem = PNSystemType(
        copy(inspiral.u[1]);
        λ₁=λ₁(p),
        λ₂=λ₂(p),
        PNOrder
    )
    n_modes = (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2
    h = Array{Complex{eltype(inspiral)}}(undef, n_modes, length(inspiral))
    h, pnsystem
end

"""
    coorbital_waveform!(storage, inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights in the co-orbital frame for the given
`inspiral` output by [`orbital_evolution`](@ref), using pre-allocated storage.

The storage is assumed to be the object returned from
[`coorbital_waveform_computation_storage`](@ref).  Other arguments are the same as in
[`coorbital_waveform`](@ref).
"""
function coorbital_waveform!(
    storage, inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(storage[2])
)
    h, pnsystem = storage
    @assert length(pnsystem.state) == length(inspiral.u[1])
    @assert length(inspiral) == size(h, 2)
    @assert (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2 == size(h, 1)
    @inbounds @fastmath for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        h!(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
    end
    h
end

"""
    coorbital_waveform(inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights in the co-orbital frame for the given
`inspiral` output by [`orbital_evolution`](@ref).

The `PNOrder` defaults to the one used to compute `inspiral`, but may be changed by passing
the keyword argument.

See also [`inertial_waveform`](@ref) for the waveform in the inertial frame.

!!! tip
    If you need this waveform at a different set of times `t′` than is currently present in
    `inspiral.t`, you should use the built-in interpolation capabilities of `inspiral`
    *first*, as in `inspiral′ = inspiral(t′)`, rather than interpolating the results of this
    function.  Or, perhaps better yet, you could select the times when calling
    `orbital_evolution` by using the `saveat` or `saves_per_orbit` keyword argument to that
    function.  These approaches will be more accurate, faster, and require less memory than
    interpolating the result of this function.

The mode weights are given starting at `ℓₘᵢₙ` (which must satisfy `0 ≤ ℓₘᵢₙ ≤ 2`) and
extending through `ℓₘₐₓ`.  The waveform is returned as a 2-dimensional `Array`, in which the
first index varies over the mode index from `(ℓ, m) = (ℓₘᵢₙ, -ℓₘᵢₙ)` to `(ℓ, m) = (ℓₘₐₓ,
ℓₘₐₓ)`, with `m` varying most rapidly, and the second index varying over the time steps.

By default, the waveform is returned in the co-orbital frame — which is somewhat like the
co-rotating frame.  In particular, the modes of non-precessing systems vary slowly, over
inspiral timescales; modes of precessing systems still vary on orbital timescales, though
even this variation could be factored out.  If `inertial=true` is passed, the waveform is
instead transformed to the inertial frame, resulting in the oscillatory behavior we usually
expect from a waveform.
"""
function coorbital_waveform(
    inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(_pnsystem(inspiral))
)
    storage = coorbital_waveform_computation_storage(
        inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder
    )
    coorbital_waveform!(
        storage, inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder
    )
end


"""
    inertial_waveform_computation_storage(inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Construct storage needed to compute waveforms in the inertial frame without allocations.

This returns the storage for the waveforms themselves, storage used for computing the Wigner
𝔇 matrices, and for "in-place" multiplication.  The returned quantity can just be passed as
the first argument to [`inertial_waveform!`](@ref) without being unpacked.

The meaning of the arguments is the same as in [`inertial_waveform`](@ref).
"""
function inertial_waveform_computation_storage(
    inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(_pnsystem(inspiral))
)
    h, pnsystem = coorbital_waveform_computation_storage(inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder)
    (D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ) = Dprep(ℓₘₐₓ, eltype(inspiral))
    hᵢ = similar(h[:,1])
    h, pnsystem, D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ, hᵢ
end

"""
    inertial_waveform!(storage, inspiral; [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights in the inertial frame for the given
`inspiral` output by [`orbital_evolution`](@ref), using pre-allocated storage.

The storage is assumed to be the object returned from
[`inertial_waveform_computation_storage`](@ref).  Other arguments are the same as in
[`inertial_waveform`](@ref).
"""
function inertial_waveform!(
    storage, inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(storage[2])
)
    h, pnsystem, D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ, hᵢ = storage
    @assert length(pnsystem.state) == length(inspiral.u[1])
    @assert length(inspiral) == size(h, 2)
    @assert (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2 == size(h, 1)
    @inbounds @fastmath for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        h!(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
        D!(D, R(pnsystem), ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
        f′ = Yiterator(hᵢ, ℓₘₐₓ, ℓₘᵢₙ, 1)
        f = Yiterator(h[:, iₜ], ℓₘₐₓ, ℓₘᵢₙ, 1)
        𝔇 = Diterator(D, ℓₘₐₓ, ℓₘᵢₙ)
        for (f′ˡ, fˡ, 𝔇ˡ) in zip(f′, f, 𝔇)
            mul!(f′ˡ, 𝔇ˡ, fˡ)
        end
        h[:, iₜ] .= hᵢ
    end
    h
end

"""
    inertial_waveform(inspiral, [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [PNOrder])

Evaluate the post-Newtonian waveform mode weights in the inertial frame for the given
`inspiral` output by [`orbital_evolution`](@ref).

The inertial frame is the one in which inertial observers are found, so this waveform is
more like one that actual observers would detect.  This function transforms the waveform
from the co-orbital frame — which is the one in which PN expressions are provided.

See [`coorbital_waveform`](@ref) for details about the other arguments.
"""
function inertial_waveform(
    inspiral; ℓₘᵢₙ=2, ℓₘₐₓ=8, PNOrder=pn_order(_pnsystem(inspiral))
)
    storage = inertial_waveform_computation_storage(
        inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder
    )
    inertial_waveform!(
        storage, inspiral; ℓₘᵢₙ, ℓₘₐₓ, PNOrder
    )
end
