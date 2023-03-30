"""
    coorbital_waveform(inspiral, [ℓₘᵢₙ=2], [ℓₘₐₓ=8], [inertial=false], [PNOrder])

Evaluate the post-Newtonian waveform mode weights for the given `inspiral` output by
[`orbital_evolution`](@ref).

See also [`inertial_waveform`](@ref) for the waveform in the inertial frame.

!!! tip

    If you need this waveform at a different set of times `t′` than is currently present in
    `inspiral.t`, you should use the built-in interpolation capabilities of `inspiral`
    *first*, as in `inspiral′ = inspiral(t′)`, rather than interpolating the results of this
    function.  Or, perhaps better yet, you could select the times when calling
    `orbital_evolution` by using the `saveat` keyword argument to that function.  These
    approaches will be more accurate, faster, and require less memory.

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
function coorbital_waveform(
    inspiral::SciMLBase.ODESolution; ℓₘᵢₙ=2, ℓₘₐₓ=8,
    inertial=false, PNOrder=pn_order(inspiral.prob.p)
)
    coorbital_waveform(
        inspiral, inspiral.prob.p; ℓₘᵢₙ, ℓₘₐₓ, inertial, PNOrder
    )
end
function coorbital_waveform(
    inspiral::SciMLBase.DiffEqArray; ℓₘᵢₙ=2, ℓₘₐₓ=8,
    inertial=false, PNOrder=pn_order(inspiral.p)
)
    coorbital_waveform(
        inspiral, inspiral.p; ℓₘᵢₙ, ℓₘₐₓ, inertial, PNOrder
    )
end
function coorbital_waveform(
    inspiral, p; ℓₘᵢₙ, ℓₘₐₓ, inertial, PNOrder
)
    @assert 0 ≤ ℓₘᵢₙ ≤ 2
    @assert ℓₘᵢₙ ≤ ℓₘₐₓ
    PNSystemType = typeof(p).name.wrapper
    pnsystem = PNSystemType(
        copy(inspiral.u[1]);
        λ₁=λ₁(p),
        λ₂=λ₂(p),
        PNOrder
    )
    n_modes = (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2
    h = Array{Complex{eltype(inspiral)}}(undef, n_modes, length(inspiral))
    if inertial
        (D, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ) = Dprep(ℓₘₐₓ, eltype(inspiral))
        hᵢ = similar(h[:,1])
    end
    @inbounds @fastmath for iₜ ∈ eachindex(inspiral)
        pnsystem.state .= inspiral.u[iₜ]
        h!(@view(h[:, iₜ]), pnsystem; ℓₘᵢₙ, ℓₘₐₓ)
        if inertial
            D!(D, R(pnsystem), ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
            #@show ℓₘᵢₙ ℓₘₐₓ size(D) size(h[:, iₜ]) size(hᵢ)
            f′ = Yiterator(hᵢ, ℓₘₐₓ, ℓₘᵢₙ, 1)
            f = Yiterator(h[:, iₜ], ℓₘₐₓ, ℓₘᵢₙ, 1)
            𝔇 = Diterator(D, ℓₘₐₓ, ℓₘᵢₙ)
            for (f′ˡ, fˡ, 𝔇ˡ) in zip(f′, f, 𝔇)
                mul!(f′ˡ, 𝔇ˡ, fˡ)
            end
            #mul!(hᵢ, D, @view(h[:, iₜ]))
            h[:, iₜ] .= hᵢ
        end
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
    inspiral::SciMLBase.ODESolution; ℓₘᵢₙ=2, ℓₘₐₓ=8,
    PNOrder=pn_order(inspiral.prob.p)
)
    coorbital_waveform(inspiral; ℓₘᵢₙ, ℓₘₐₓ, inertial=true, PNOrder)
end
function inertial_waveform(
    inspiral::SciMLBase.DiffEqArray; ℓₘᵢₙ=2, ℓₘₐₓ=8,
    PNOrder=pn_order(inspiral.p)
)
    coorbital_waveform(inspiral; ℓₘᵢₙ, ℓₘₐₓ, inertial=true, PNOrder)
end
