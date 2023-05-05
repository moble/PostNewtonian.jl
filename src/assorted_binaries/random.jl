@doc raw"""
    rand([rng=GLOBAL_RNG], pnclass; [v=0.01], [PNOrder=typemax(Int)])

Generate a random PNSystem, specifically of class `pnclass` (one of `BBH`, `BHNS`, or
`NSNS`), with parameters chosen to be consistent with typical LIGO/Virgo/KAGRA searches as
of this writing.

# Examples
```julia-repl
julia> pnsystem = rand(BBH);

julia> inspiral = orbital_evolution(pnsystem);

```

# Parameter ranges

Choosing the space of parameters from which to randomly select values is nontrivial.  The
most relevant numbers are of two types from LIGO/Virgo/KAGRA: (1) the range of parameters
over which systems are *searched for*, and (2) the range of parameters used as priors in
parameter estimation.

Appendix B1 of LIGO's [GWTC-1 catalog paper](https://arxiv.org/abs/1811.12907) says that the
"spin vectors are assumed to be isotropic on the sphere and uniform in spin magnitude", with
two choices of uniform magnitude: ``a_i \leq 0.89`` and ``a_i \leq 0.05``.  The
dimensionless tidal deformabilities ``\Lambda_i`` of each NS are assumed to be jointly
uniform within ``0 \leq \Lambda_i \leq 5000``.

The more current [GWTC-3 paper](https://arxiv.org/abs/2111.03606) searches with the
following:

> The PyCBC-BBH analysis focuses on a region ranging in primary component mass from 5M to
> 350M, with mass ratios from 1/3 to 1, and effective spins ranging from ``\chi_\mathrm{eff}
> = −0.998`` to ``0.998``.

That paper's parameter estimation uses uniform priors over spin magnitudes and redshifted
component masses; isotropic spin orientation, sky location and binary orientation;
mass-ratio ``q \in [0.05, 1]`` (sometimes extending down to ``q=0.02``).  It seems that the
actual ranges of the spin magnitudes and masses are restricted based on initial guesses,
which may be expanded to ensure that the posteriors lie entirely within the priors.  The
distance priors are uniform in ``D_\mathrm{L}^2``, with some adjustments due to cosmology.

For the purposes of this package, we're not too interested in absolute scales — specifically
distance or total mass — so we'll just choose the total mass to be ``1`` at the initial
frequency, then sample based on mass ratio and spin magnitudes.

Given these, it seems like the current state of the art would be well covered by choosing a
random `q` ∈ [0.05, 1], `χᵢ` ∈ [0, 0.998], and `Λ̂ᵢ` ∈ [0, 5000], with isotropic choices for
orientations.
"""
Base.rand(pnclass::Type{P}; v::T=0.2, PNOrder=typemax(Int)) where {P<:PNSystem, T} =
    rand(GLOBAL_RNG, pnclass; v, PNOrder)

function Base.rand(
    rng::AbstractRNG, pnclass::Type{P}; v::T=0.2, PNOrder=typemax(Int)
) where {P<:PNSystem, T}
    qₘᵢₙ = T(big"0.05")  # Note that we're using q≤1 here for consistency with LIGO
    χₘₐₓ = T(big"0.998")
    Λ̂ₘₐₓ = T(big"5000.0")
    q = rand(rng, T) * (1-qₘᵢₙ) + qₘᵢₙ
    M₁ = 1/(q+1)
    M₂ = q/(q+1)
    χ⃗₁ = χₘₐₓ * rand(rng, T) * normalize(randn(rng, QuatVec{T}))
    χ⃗₂ = χₘₐₓ * rand(rng, T) * normalize(randn(rng, QuatVec{T}))
    R = randn(rng, Rotor{T})
    Λ̂₁ = Λ̂ₘₐₓ * rand(rng, T)
    Λ̂₂ = Λ̂ₘₐₓ * rand(rng, T)
    pnclass(;M₁, M₂, χ⃗₁, χ⃗₂, R, v, Λ̂₁, Λ̂₂, PNOrder)
end
