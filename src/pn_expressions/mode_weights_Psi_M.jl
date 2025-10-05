@doc raw"""
    Ψ_M!(Ψ_M, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))
    mode_weights_Ψ_M!(Ψ_M, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))

Compute mode weights of Moreschi supermomentum Ψ_M emitted by `pn` system, modifying `Ψ_M` in place.

!!! note
    This is a low-level function; you probably don't want to use this directly.  See
    [`coorbital_waveform`](@ref) or [`inertial_waveform`](@ref) for more user-friendly
    functions.

These modes are computed in the "co-orbital" frame, in which the larger object lies on the
positive ``x`` axis, the smaller lies on the negative ``x`` axis, and the instantaneous
angular velocity is in the positive ``z`` direction.

The modes are stored in `Ψ_M` in order of increasing ``ℓ`` and increasing ``m``, with ``m``
iterating fastest, all the way up to the highest available mode, ``(8,8)``.

The Moreschi supermomentum has spin weight 0, and the ``(ℓ,m)=(0,0)`` mode is already nonzero.
Nonetheless a user may request some `ℓₘᵢₙ` to the smallest ``ℓ`` value that
should be present in the output data.

These results come from Appendix A of [Mitman, Stein, Boyle, et al.
(2022)](https://arxiv.org/abs/2208.04356).
"""
@pn_expression 2 function Ψ_M!(Ψ_M, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))
    Ψ_M .= false  # Set everything to 0 just to be safe

    # Eq. (46). Recall that x=(v/c)^2
    # Note that Eq. (46) shouldn't have the 1/R^2 in its first term.
    Ψ_M₀ = 2M*ν * (v/c)^2 * √(π/4)

    # ell=0
    if ℓₘₐₓ ≥ 0
        # Eq. (46) (the leading -M) and Eq. (A1a) of MSB+ (2022).
        Ψ_M[Yindex(0,0,ℓₘᵢₙ)] = -M + Ψ_M₀ * @pn_expansion(
            1
            + (v/c)^2 * (-3//4-ν/12)
            + (v/c)^4 * (-27//8+19ν/8-ν^2/24)
            + (v/c)^6 * (-675//64+(34445//576-205π^2/96)ν-155ν^2/96-35ν^3/5184)
        )
    end

    # ell=1
    if ℓₘₐₓ ≥ 1
        Ψ_M[Yindex(1,1,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1c) of MSB+ (2022)
            (v/c)^6 * (δ * ν * 464//35 * √(2//3))
        )
    end

    # ell=2
    if ℓₘₐₓ ≥ 2
        Ψ_M[Yindex(2,0,ℓₘᵢₙ)] = Ψ_M₀ * (2√5)/7 * @pn_expansion(
            # Eq. (A1g) of MSB+ (2022)
            1
            + (v/c)^2 * (-4075/4032 + 67ν/48)
            + (v/c)^4 * (-151877213//67060224 - 123815ν/44352 + 205ν^2/352)
            + (v/c)^5 * (-253/336 + 253ν/84)π
            + (v/c)^6 * (-4397711103307//532580106240
                         + (700464542023//13948526592 - 205π^2/96)ν
                         + 69527951ν^2/166053888
                         + 1321981ν^3/5930496)
        )
    end

    # ell=3
    if ℓₘₐₓ ≥ 3
        Ψ_M[Yindex(3,1,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1k) of MSB+ (2022)
            (v/c)^6 * (484δ*ν/(15√21))
        )
        Ψ_M[Yindex(3,3,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1i) of MSB+ (2022)
            (v/c)^6 * (-44δ*ν/(27√35))
        )
    end

    # ell=4
    if ℓₘₐₓ ≥ 4
        Ψ_M[Yindex(4,0,ℓₘᵢₙ)] = Ψ_M₀ * (1//42) * @pn_expansion(
            # Eq. (A1q) of MSB+ (2022)
            1
            + (v/c)^2 * (-180101//29568 + 27227ν/1056)
            + (v/c)^4 * (2201411267//158505984
                         - 34829479ν/432432 + 844951ν^2/27456)
            + (v/c)^5 * (-13565/1232 + 13565ν/308)π
            + (v/c)^6 * (15240463356751//781117489152
                         - (1029744557245//27897053184 + 205π^2/96)*ν
                         - 4174614175ν^2/36900864
                         + 221405645ν^3/11860992)
        )
        Ψ_M[Yindex(4,4,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1m) of MSB+ (2022)
            (v/c)^5 * (-4𝒾//3 * √(2//35) * ν)
        )
    end

    # ell=5
    if ℓₘₐₓ ≥ 5
        Ψ_M[Yindex(5,1,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1w) of MSB+ (2022)
            (v/c)^6 * (52//21*√(2//165)*δ*ν)
        )
        Ψ_M[Yindex(5,3,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1u) of MSB+ (2022)
            (v/c)^6 * (4δ*ν/(27√385))
        )
        Ψ_M[Yindex(5,5,ℓₘᵢₙ)] = Ψ_M₀ * @pn_expansion(
            # Eq. (A1s) of MSB+ (2022)
            (v/c)^6 * (-36δ*ν/(5√77))
        )
    end

    # ell=6
    if ℓₘₐₓ ≥ 6
        Ψ_M[Yindex(6,0,ℓₘᵢₙ)] = Ψ_M₀ * (-4195/(177408√13)) * @pn_expansion(
            # Eq. (A1y) of MSB+ (2022)
            + (v/c)^2 * (1 - 3612ν/839)
            + (v/c)^4 * (-45661561//6342840 + 101414ν/2517 - 48118ν^2/839)
            + (v/c)^5 * (1248//839 - 4992ν/839)π
            + (v/c)^6 * (3012132889099//144921208320
                         - 27653500031ν/191694720
                         + 1317967427ν^2/4107744
                         - 24793657ν^3/342312)
        )
    end

    # ell=7
    # Nothing so far

    # ell=8
    if ℓₘₐₓ ≥ 8
        Ψ_M[Yindex(8,0,ℓₘᵢₙ)] = Ψ_M₀ * (75601/(8895744√17)) * @pn_expansion(
            # Eq. (A1aa) of MSB+ (2022)
            + (v/c)^4 * (1 - 452070ν/75601 + 733320ν^2/75601)
            + (v/c)^6 * (-265361599//33869248
                         + 18177898147ν/321757856
                         - 722521125ν^2/5745676
                         + 261283995ν^3/2872838)
        )
    end

    # Spin terms
    if ℓₘₐₓ ≥ 0
        Ψ_M[Yindex(0,0,ℓₘᵢₙ)] += Ψ_M₀ * @pn_expansion(
            # Eq. (A1b) of MSB+ (2022)
            + (v/c)^3 * (14//3 * Sₗ + 2δ * Σₗ) / M^2
            + (v/c)^4 * (-(16S⃗·S⃗ + 3Σ⃗·Σ⃗ + 32Sₗ^2 + 9Σₗ^2)/12
                         - 4//3 * δ * (S⃗·Σ⃗ + 2Sₗ*Σₗ)
                         + 4//3 * ν * (Σ⃗·Σ⃗ + 2Σₗ^2)) / M^4
        )
    end
    if ℓₘₐₓ ≥ 2
        Ψ_M[Yindex(2,0,ℓₘᵢₙ)] += Ψ_M₀ * 2//7 * √5 * @pn_expansion(
            # Eq. (A1h) of MSB+ (2022)
            + (v/c)^3 * (16//3 * Sₗ + 419//160 * δ*Σₗ) / M^2
            + (v/c)^4 * (-(128S⃗·S⃗ + 24Σ⃗·Σ⃗ + 256Sₗ^2 + 75Σₗ^2) / 96
                         - 4//3 * δ * (S⃗·Σ⃗ + 2Sₗ*Σₗ)
                         + 4//3 * ν * (Σ⃗·Σ⃗ + 2Σₗ^2)) / M^4
        )
        Ψ_M[Yindex(2,1,ℓₘᵢₙ)] += Ψ_M₀ * 61/(14√30) * @pn_expansion(
            # Eq. (A1f) of MSB+ (2022)
            + (v/c)^3 * (- (Sₙ-𝒾*Sλ) - 375//488 * δ * (Σₙ-𝒾*Σλ)) / M^2
            + (v/c)^4 * (10(3Sₗ*(Sₙ-𝒾*Sλ)+Σₗ*(Σₙ-𝒾*Σλ))
                         +15δ*((Sₙ-𝒾*Sλ)*Σₗ+Sₗ*(Σₙ-𝒾*Σλ))
                         -30ν*Σₗ*(Σₙ-𝒾*Σλ))/61M^4
        )
    end
    # ell=3
    # nothing so far
    if ℓₘₐₓ ≥ 4
        Ψ_M[Yindex(4,0,ℓₘᵢₙ)] += Ψ_M₀ * 1//42 * @pn_expansion(
            # Eq. (A1r) of MSB+ (2022)
            + (v/c)^3 * (10Sₗ + 57//8 * δ * Σₗ) / M^2
            + (v/c)^4 * (-(64S⃗·S⃗ + 12Σ⃗·Σ⃗ + 128Sₗ^2 + 41Σₗ^2) / 48
                         - 4//3 * δ * (S⃗·Σ⃗ + 2Sₗ*Σₗ)
                         + 4//3 * ν * (Σ⃗·Σ⃗ + 2Σₗ^2)) / M^4
        )
        Ψ_M[Yindex(4,1,ℓₘᵢₙ)] += Ψ_M₀ * 13/(56√5) * @pn_expansion(
            # Eq. (A1p) of MSB+ (2022)
            + (v/c)^3 * (- (Sₙ-𝒾*Sλ) - 34//39 * δ * (Σₙ-𝒾*Σλ)) / M^2
            + (v/c)^4 * (10(3Sₗ*(Sₙ-𝒾*Sλ)+Σₗ*(Σₙ-𝒾*Σλ))/3
                         +5δ*((Sₙ-𝒾*Sλ)*Σₗ+Sₗ*(Σₙ-𝒾*Σλ))
                         -10ν*Σₗ*(Σₙ-𝒾*Σλ))/39M^4
        )
    end

    # Impose reality of Ψ_M
    for ℓ in 0:ℓₘₐₓ
        for m in 1:ℓ
            Ψ_M[Yindex(ℓ,-m,ℓₘᵢₙ)] = ifelse(isodd(m), -1, 1) * conj(Ψ_M[Yindex(ℓ,m,ℓₘᵢₙ)])
        end
    end

    Ψ_M
end
const mode_weights_Ψ_M! = Ψ_M!
