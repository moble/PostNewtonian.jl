"""
    waveform(inspiral, [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [inertial=false], [PNOrder])

Evaluate the post-Newtonian waveform mode weights for the given `inspiral`.

!!! note

    If you need this waveform at a different set of times `t′` than is currently present in
    `inspiral.t`, you should use the built-in interpolation capabilities of `inspiral`
    *first*, as in `inspiral′ = inspiral(t′)`, rather than interpolating the results of this
    function.

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

The `PNOrder` defaults to the one used to compute `inspiral`, but may be given differently
here.
"""
function waveform(
    inspiral::SciMLBase.ODESolution; ℓₘᵢₙ=2, ℓₘₐₓ=8,
    inertial=false, PNOrder=pn_order(inspiral.prob.p)
)
    @assert 0 ≤ ℓₘᵢₙ ≤ 2
    @assert ℓₘᵢₙ ≤ ℓₘₐₓ
    PNSystemType = typeof(inspiral.prob.p).name.wrapper
    pnsystem = PNSystemType(
        copy(inspiral.u[1]);
        λ₁=λ₁(inspiral.prob.p),
        λ₂=λ₂(inspiral.prob.p),
        PNOrder
    )
    n_modes = (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2
    h = Array{Complex{eltype(inspiral)}}(undef, n_modes, length(inspiral))
    @inbounds for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        h!(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
    end
    h
end
