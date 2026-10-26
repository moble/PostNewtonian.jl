@doc raw"""
    h!(h, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))
    mode_weights!(h, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))

Compute mode weights of gravitational waves emitted by `pn` system, modifying `h` in place.

!!! note
    This is a low-level function; you probably don't want to use this directly.  See
    [`coorbital_waveform`](@ref) or [`inertial_waveform`](@ref) for more user-friendly
    functions.

These modes are computed in the "co-orbital" frame, in which the larger object lies on the
positive ``x`` axis, the smaller lies on the negative ``x`` axis, and the instantaneous
angular velocity is in the positive ``z`` direction.

The modes are stored in `h` in order of increasing ``ℓ`` and increasing ``m``, with ``m``
iterating fastest, all the way up to the highest available mode, ``(8,8)``.

Because gravitational waves have spin weight -2, the ``(ℓ,m)=(0,0)``, ``(1,-1)``, ``(1,0)``,
and ``(1,1)`` modes are always 0.  By default, we assume that these modes are nonetheless
included in `h`.  If that is not the case, set `ℓₘᵢₙ` to the smallest ``ℓ`` value that
should be present in the output data — `ℓₘᵢₙ=2` being the most reasonable alternative.

These results come most directly from Eqs. (A5) of [Boyle et al.
(2014)](https://arxiv.org/abs/1409.4431), with the exception of errors in the 2PN spin-spin
terms, in which cases we must multiply by ``\nu/2`` and make the substitutions $S_1 \mapsto
S_0^+$ and $S_2 \mapsto S_0^-$.  In turn, those expressions are synthesized from the
following:  Non-spinning terms are taken from [Blanchet
(2014)](https://doi.org/10.12942/lrr-2014-2), except for the highest-pN term in the (2,±2)
mode, which are taken from [Blanchet et al.  (2023)](https://arxiv.org/abs/2304.11186), and
the ``m=0`` modes, which are taken from [Favata (2008)](https://arxiv.org/abs/0812.0069).
The 1PN spin-orbit term is from Eq. (3.22d) of [Kidder
(1995)](https://link.aps.org/doi/10.1103/PhysRevD.52.821).  The 1.5PN spin-orbit term is
from Eq. (3.22f) of Kidder (1995) and Eq. (F15b) of [Will and Wiseman
(1996)](https://link.aps.org/doi/10.1103/PhysRevD.54.4813).  The 2PN spin-orbit term is from
Eq. (4.13) of [Buonanno, Faye, Hinderer
(2013)](https://link.aps.org/doi/10.1103/PhysRevD.87.044009), while the 2PN spin-spin term
is from Eq. (4.15) of that reference.
"""
@pn_expression 2 function h!(h, pnsystem; ℓₘᵢₙ=0, ℓₘₐₓ=typemax(Int))
    h .= false  # Set everything to 0 just to be safe

    h₀ = 2ν * (v/c)^2 * √(16π/5)

    # ell=2
    if ℓₘₐₓ ≥ 2
        h[Yindex(2,0,ℓₘᵢₙ)] = h₀ * (-5/(14√6)) * @pn_expansion(
            # Eq. (4.3a) of Favata (2008)
            1
            + (v/c)^2 * (-4075//4032 + 67ν/48)
            + (v/c)^4 * (-151877213//67060224 - 123815ν/44352 + 205ν^2/352)
            + (v/c)^5 * (-253//336 + 253ν/84)π
            + (v/c)^6 * (
                -4397711103307//532580106240
                + (700464542023//13948526592 - 205π^2/96)ν
                + 69527951ν^2/166053888
                + 1321981ν^3/5930496
            )
        )
        h[Yindex(2,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^1 * (𝒾 * δ / 3)
            + (v/c)^3 * (𝒾 * δ * (-17 + 20ν) / 84)
            + (v/c)^4 * ((δ * (1 + 2𝒾*π + 4ln2)) / 6)
            + (v/c)^5 * (𝒾 * δ * (-172 + ν * (-2036 + 237ν)) / 1512)
            + (v/c)^6 * (δ*(-34𝒾*π - 17*(1 + 4ln2) + 2ν * (353 + 6𝒾*π + 12ln2)) / 168)
        )
        h[Yindex(2,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            1
            + (v/c)^2 * ((-107 + 55ν)/42)
            + (v/c)^3 * (2π)
            + (v/c)^4 * ((-2173 + ν * (-7483 + 2047ν)) / 1512)
            + (v/c)^5 * (-24𝒾*ν + ((-107 + 34ν)π) / 21)
            + (v/c)^6 * (
                (27027409//646800) - 856γₑ/105 + (ν*(-834555 + ν*(-729396 + 114635ν))) / 99792
                + 41ν * π^2 / 96 + (2π*(214𝒾 + 35π))/105 - 1712ln2/105
                - (856//105)*ln(v)
            )
            + (v/c)^7 * ((-2𝒾 * ν * (-501655 + 24396ν) + 15*(-2173 + 2ν*(-2495 + 560ν))π) / 11340)
            # Eq. (6.17) of Blanchet et al. (2023)
            + (v/c)^8 * (
                - 846557506853//12713500800 + 45796γₑ/2205 - 22898𝒾*π/2205 - 107π^2/63 + 45796(2ln2+ln(v))/2205
                + (-336005827477//4237833600 + 15284γₑ/441 - 219314𝒾*π/2205 - 9755*π^2/32256 + 15284(2ln2+ln(v))/441)ν
                + (256450291//7413120 - 1025*π^2/1008)ν^2 - 81579187ν^3/15567552 + 26251249ν^4/31135104
            )
        )
    end

    # ell=3
    if ℓₘₐₓ ≥ 3
        h[Yindex(3,0,ℓₘᵢₙ)] = h₀ * @pn_expansion((v/c)^5 * (-2𝒾 * √(6//7) * ν / 5))
        h[Yindex(3,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^1 * (𝒾 * δ / 12√14)
            + (v/c)^3 * (-𝒾 * δ * (4 + ν) / 18√14)
            + (v/c)^4 * (δ * (7 + 5𝒾*π + 10ln2) / 60√14)
            + (v/c)^5 * (-𝒾 * δ * (-607 + ν*(272 + 247ν)) / 2376√14)
            + (v/c)^6 * (
                δ * (-5𝒾 * (16 + 7ν)π + 2*(-56 + ν - 5*(16 + 7ν)*ln2)) / 360√14
            )
            + (v/c)^7 * (
                𝒾 * δ / 12√14 * (
                    (10753397//1513512) - 2ln2*((212//105) + ln2) - (26//21)*γₑ + (π^2/6)
                    - 2𝒾*π*((41//105) + ln2) + (ν/8)*(-(1738843//19305) + (41//8)*π^2)
                    + (327059//30888)*ν^2 - (17525//15444)*ν^3
                    + ln(v) * (-26//21)
                )
            )
        )
        h[Yindex(3,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^2 * (√(5//7) * (1 - 3ν) / 3)
            + (v/c)^4 * ((-193 + (145 - 73ν)*5ν) / 54√35)
            + (v/c)^5 * ((-15𝒾 + 66𝒾*ν + 10π - 30π*ν) / 3√35)
            + (v/c)^6 * ((-1451 + (-17387 + 3*(33342 - 5341ν)ν)ν) / 2376√35)
        )
        h[Yindex(3,3,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^1 * (-3𝒾 * √(15//224) * δ)
            + (v/c)^3 * (-3𝒾 * √(15//56) * δ * (-2 + ν))
            + (v/c)^4 * (√(243//70) * δ * (-7 - 5𝒾*π + 10ln³╱₂) / 4)
            + (v/c)^5 * (-𝒾 * √(3//542080) * δ * (369 + (-3676 + 887ν)ν))
            + (v/c)^6 * (δ * (-3645𝒾 * (-8 + 3ν)π + (40824 - 96206ν + 7290*(-8 + 3ν)ln³╱₂)) / 216√210)
            + (v/c)^7 * (
                ((-3𝒾)/4) * √(15//14) * δ * (
                    (19388147//280280) + (492//35)ln³╱₂ - 18*ln³╱₂^2 - (78//7)γₑ + (3//2)π^2
                    + 6𝒾 * π * (-41//35 + 3ln³╱₂)
                    + (-7055//3432 + 41//64 * π^2)ν - (318841//17160)ν^2 + (8237//2860)ν^3
                    + ln(v) * (-(39//7) * 4ln2)
                )
            )
        )
    end

    # ell=4
    if ℓₘₐₓ ≥ 4
        h[Yindex(4,0,ℓₘᵢₙ)] = h₀ * (-1 / 504√2) * @pn_expansion(
            # Eq. (4.3b) of Favata (2008)
            1
            + (v/c)^2 * (-180101//29568 + 27227ν/1056)
            + (v/c)^4 * (2201411267//158505984 - 34829479ν/432432 + 844951ν^2/27456)
            + (v/c)^5 * (-13565//1232 + 13565ν/308)π
            + (v/c)^6 * (
                15240463356751//781117489152
                + (-1029744557245//27897053184 - 205π^2/96)ν
                - 4174614175ν^2/36900864
                + 221405645ν^3/11860992
            )
        )
        h[Yindex(4,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^3 * (𝒾 * δ * (1 - 2ν) / 84√10)
            + (v/c)^5 * (-𝒾 * δ * (404 + (-1011 + 332ν)ν) / 11088√10)
            + (v/c)^6 * (δ * (64 - 1661ν - 30𝒾*(-1 + 2ν)π + 60*(1 - 2ν)ln2) / 2520√10)
        )
        h[Yindex(4,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^2 * (√5 * (1 - 3ν) / 63)
            + (v/c)^4 * ((-1311 + 5*(805 - 57ν)ν) / 4158√5)
            + (v/c)^5 * ((-21𝒾 + ν*(84𝒾 - 30π) + 10π) / 63√5)
            + (v/c)^6 * ((9342351 + 7ν*(-5460759 + 115ν*(34822 + 3363ν))) / 22702680√5)
        )
        h[Yindex(4,3,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^3 * (9𝒾 * δ * (-1 + 2ν) / 4√70)
            + (v/c)^5 * (3𝒾 * δ * (468 + (-1267 + 524ν)ν) / 176√70)
            + (v/c)^6 * (δ * (-5184 + 16301ν + 2430𝒾*(-1 + 2ν)π + 4860*(1 - 2ν)ln³╱₂) / 360√70)
        )
        h[Yindex(4,4,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^2 * (8 * √(5//7) * (-1 + 3ν) / 9)
            + (v/c)^4 * (4 * (1779 + 5ν*(-1273 + 525ν)) / 297√35)
            + (v/c)^5 * ((160*(-1 + 3ν)π + 𝒾*(336 - 1193ν + 320*(-1 + 3ν)ln2)) / 9√35)
            + (v/c)^6 * ((-9618039 + 7ν*(9793071 + 5ν*(-3231338 + 678291ν))) / 405405√35)
        )
    end

    # ell=5
    if ℓₘₐₓ ≥ 5
        h[Yindex(5,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^3 * (𝒾 * δ * (1 - 2ν) / 288√385)
            + (v/c)^5 * (-𝒾 * δ * (179 + 4*(-88 + ν)ν) / 11232√385)
            + (v/c)^6 * (δ * (181 - 70𝒾*(-1 + 2ν)π + 140ln2 - 28ν*(313 + 10ln2)) / 20160√385)
        )
        h[Yindex(5,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^4 * ((2 + 10*(-1 + ν)ν) / 27√55)
            + (v/c)^6 * ((-3911 + 7ν*(3079 + 35ν*(-118 + 33ν))) / 12285√55)
        )
        h[Yindex(5,3,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^3 * (9𝒾 * √(3//110) * δ * (-1 + 2ν) / 32)
            + (v/c)^5 * (3𝒾 * √(3//110) * δ * (207 + 8ν*(-58 + 11ν)) / 416)
            + (v/c)^6 * (δ * (-395847 + 1171828ν + 153090𝒾*(-1 + 2ν)π - 306180*(-1 + 2ν)*ln³╱₂) / 60480√330)
        )
        h[Yindex(5,4,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^4 * ((-32 - 160*(-1 + ν)ν) / 9√165)
            + (v/c)^6 * (16*(4451 - 7ν*(3619 + 5ν*(-1042 + 339ν))) / 4095√165)
        )
        h[Yindex(5,5,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^3 * (-625𝒾 * δ * (-1 + 2ν) / 96√66)
            + (v/c)^5 * (-625𝒾 * δ * (263 + 16ν*(-43 + 16ν)) / 3744√66)
            + (v/c)^6 * (δ * (565625 - 1481676ν - 218750𝒾*(-1 + 2ν)π + 437500*(-1 + 2ν)ln⁵╱₂) / 6720√66)
        )
    end

    # ell=6
    if ℓₘₐₓ ≥ 6
        h[Yindex(6,0,ℓₘᵢₙ)] = h₀ * (4195/(1419264√273)) * @pn_expansion(
            # Eq. (4.3c) of Favata (2008)
            + (v/c)^2 * (1 - 3612ν/839)
            + (v/c)^4 * (-45661561//6342840 + 101414ν/2517 - 48118ν^2/839)
            + (v/c)^5 * (1248//839 - 4992ν/839)π
            + (v/c)^6 * (
                3012132889099//144921208320
                - 27653500031ν/191694720
                + 1317967427ν^2/4107744
                - 24793657ν^3/342312
            )
        )
        h[Yindex(6,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (𝒾 * δ * (-1 + ν) * (-1 + 3 * ν) / 8316√26)
        )
        h[Yindex(6,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^4 * ((2 + (-1 + ν) * 10ν) / 297√65)
            + (v/c)^6 * ((-81 + (59 + (-64 + 7ν)ν) * 7ν) / 2079√65)
        )
        h[Yindex(6,3,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (-81𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 616√65)
        )
        h[Yindex(6,4,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^4 * (-128 * √(2//39) * (1 + 5 * (-1 + ν)ν) / 495)
            + (v/c)^6 * (-64 * √(2//39) * (-93 + 7 * (71 + (-88 + 19ν)ν)ν) / 3465)
        )
        h[Yindex(6,5,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (3125𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 504√429)
        )
        h[Yindex(6,6,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^4 * (54 * (1 + 5 * (-1 + ν)ν) / 5√143)
            + (v/c)^6 * (27 * (-113 + 7 * (91 + (-128 + 39ν)ν)ν) / 35√143)
        )
    end

    # ell=7
    if ℓₘₐₓ ≥ 7
        h[Yindex(7,1,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 864864√2)
        )
        h[Yindex(7,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * ((1 - (-1 + ν)^2 * 7ν) / 3003√3)
        )
        h[Yindex(7,3,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (-243𝒾 * √(3//2) * δ * (-1 + ν) * (-1 + 3ν) / 160160)
        )
        h[Yindex(7,4,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * (128√(2//33) * (-1 + (-1 + ν)^2 * 7ν) / 1365)
        )
        h[Yindex(7,5,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (15625𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 26208√66)
        )
        h[Yindex(7,6,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * (-81√(3//143) * (-1 + (-1 + ν)^2 * 7ν) / 35)
        )
        h[Yindex(7,7,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^5 * (-16807𝒾 * √(7//858) * δ * (-1 + ν) * (-1 + 3ν) / 1440)
        )
    end

    # ell=8
    if ℓₘₐₓ ≥ 8
        h[Yindex(8,0,ℓₘᵢₙ)] = h₀ * (-75601/(213497856√119)) * @pn_expansion(
            # Eq. (4.3d) of Favata (2008)
            + (v/c)^4 * (1 - 452070ν/75601 + 733320ν^2/75601)
            + (v/c)^6 * (
                - 265361599//33869248
                + 18177898147ν/321757856
                - 722521125ν^2/5745676
                + 261283995ν^3/2872838
            )
        )
        h[Yindex(8,2,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * (-(-1 + (-1 + ν)^2 * 7ν) / 9009√85)
        )
        h[Yindex(8,4,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * ((128√(2//187) * (-1 + (-1 + ν)^2 * 7ν)) / 4095)
        )
        h[Yindex(8,6,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * ((-243√(3//17017) * (-1 + (-1 + ν)^2 * 7ν)) / 35)
        )
        h[Yindex(8,8,ℓₘᵢₙ)] = h₀ * @pn_expansion(
            (v/c)^6 * ((16384√(2//85085) * (-1 + (-1 + ν)^2 * 7ν)) / 63)
        )
    end

    # Symmetric spin terms
    if ℓₘₐₓ ≥ 2
        h[Yindex(2,0,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^4 * ((S₀⁺ₙ * S₀⁻ₙ - S₀⁺λ * S₀⁻λ) / (√6 * M^4))
        )
        h[Yindex(2,1,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^2 * (𝒾 * Σₗ / 2M^2)
            + (v/c)^4 * (𝒾 * (-86*Sₗ*δ + Σₗ*(139ν - 79)) / 42M^2)
        )
        h[Yindex(2,2,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^3 * (-(6Sₗ + 2Σₗ*δ) / 3M^2)
            + (v/c)^4 * (
                (12S₀⁺ₗ*S₀⁻ₗ + 10S₀⁺λ*S₀⁻λ - 22S₀⁺ₙ*S₀⁻ₙ - 15𝒾 * (S₀⁺ₙ*S₀⁻λ + S₀⁺λ*S₀⁻ₙ))
                / 12M^4
            )
        )
    end
    if ℓₘₐₓ ≥ 3
        h[Yindex(3,1,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^4 * (√14𝒾 * (Sₗ*δ - 5Σₗ*(3ν - 1)) / 336M^2)
        )
        h[Yindex(3,2,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^3 * (2√35 * (Sₗ + Σₗ*δ) / 21M^2)
        )
        h[Yindex(3,3,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^4 * (3√210𝒾 * (7Sₗ*δ - 3Σₗ*(3ν - 1)) / 112M^2)
        )
    end
    if ℓₘₐₓ ≥ 4
        h[Yindex(4,1,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^4 * (√10𝒾 * (Sₗ*δ - 3Σₗ*ν + Σₗ) / 336M^2)
        )
        h[Yindex(4,3,ℓₘᵢₙ)] += h₀ * @pn_expansion(
            (v/c)^4 * (9√70𝒾 * (-Sₗ*δ + 3Σₗ*ν - Σₗ) / 112M^2)
        )
    end

    # Symmetrize everything
    for ℓ in 2:ℓₘₐₓ
        for m in 1:ℓ
            h[Yindex(ℓ,-m,ℓₘᵢₙ)] = ifelse(isodd(ℓ), -1, 1) * conj(h[Yindex(ℓ,m,ℓₘᵢₙ)])
        end
    end

    # Anti-symmetric spin terms
    if ℓₘₐₓ ≥ 2
        h̃₂₀ = h₀ * @pn_expansion(
            (v/c)^2 * (√6𝒾 * Σₙ / 6M^2)
            + (v/c)^4 * (√6𝒾 * (255Sₙ*δ - Σₙ*(506ν - 45)) / 126M^2)
        )
        h[Yindex(2,0,ℓₘᵢₙ)] += h̃₂₀
        h̃₂₁ = h₀ * @pn_expansion(
            (v/c)^3 * ((4𝒾*Sλ + 25*Sₙ + 4𝒾*Σλ*δ + 13*Σₙ*δ) / 6M^2)
            + (v/c)^4 * -3 * (S₀⁺ₙ*S₀⁻ₗ + S₀⁺ₗ*S₀⁻ₙ) / (4M^4)
        )
        h[Yindex(2,1,ℓₘᵢₙ)] += h̃₂₁
        h[Yindex(2,-1,ℓₘᵢₙ)] += -conj(h̃₂₁)
        h̃₂₂ = h₀ * @pn_expansion(
            (v/c)^2 * (-(Σλ + 𝒾*Σₙ) / 2M^2)
            + (v/c)^4 * ((19*Sλ*δ + 182𝒾*Sₙ*δ - 43*Σλ*ν + 5*Σλ - 280𝒾*Σₙ*ν + 98𝒾*Σₙ) / 84M^2)
        )
        h[Yindex(2,2,ℓₘᵢₙ)] += h̃₂₂
        h[Yindex(2,-2,ℓₘᵢₙ)] += -conj(h̃₂₂)
    end
    if ℓₘₐₓ ≥ 3
        h̃₃₀ = h₀ * @pn_expansion(
            (v/c)^4 * (√42 * (-17Sλ*δ + Σλ*(35ν - 9)) / 168M^2)
        )
        h[Yindex(3,0,ℓₘᵢₙ)] += h̃₃₀
        h̃₃₁ = h₀ * @pn_expansion(
            (v/c)^3 * (√14 * (𝒾*Sλ + Sₙ + δ*(𝒾*Σλ + Σₙ)) / 21M^2)
        )
        h[Yindex(3,1,ℓₘᵢₙ)] += h̃₃₁
        h[Yindex(3,-1,ℓₘᵢₙ)] += conj(h̃₃₁)
        h̃₃₂ = h₀ * @pn_expansion(
            (v/c)^4 * (√35 * (-Σλ*(83ν - 17) + 4𝒾*Σₙ*(55ν - 13) + 25δ * (Sλ - 4𝒾*Sₙ)) / 168M^2)
        )
        h[Yindex(3,2,ℓₘᵢₙ)] += h̃₃₂
        h[Yindex(3,-2,ℓₘᵢₙ)] += conj(h̃₃₂)
        h̃₃₃ = h₀ * @pn_expansion(
            (v/c)^3 * (√210𝒾 * (Sλ + 𝒾*Sₙ + δ*(Σλ + 𝒾*Σₙ)) / 21M^2)
        )
        h[Yindex(3,3,ℓₘᵢₙ)] += h̃₃₃
        h[Yindex(3,-3,ℓₘᵢₙ)] += conj(h̃₃₃)
    end
    if ℓₘₐₓ ≥ 4
        h̃₄₀ = h₀ * @pn_expansion(
            (v/c)^4 * (√2𝒾 * (Sₙ*δ - 3Σₙ*ν + Σₙ) / 168M^2)
        )
        h[Yindex(4,0,ℓₘᵢₙ)] += h̃₄₀
        h̃₄₂ = h₀ * @pn_expansion(
            (v/c)^4 * (√5 * (-13Σλ*(3ν - 1) + 14𝒾*Σₙ*(3ν - 1) + δ*(13Sλ - 14𝒾*Sₙ)) / 168M^2)
        )
        h[Yindex(4,2,ℓₘᵢₙ)] += h̃₄₂
        h[Yindex(4,-2,ℓₘᵢₙ)] += -conj(h̃₄₂)
        h̃₄₄ = h₀ * @pn_expansion(
            (v/c)^4 * (9√35 * (-3Σλ*ν + Σλ - 𝒾*Σₙ*(3ν - 1) + δ*(Sλ + 𝒾*Sₙ)) / 56M^2)
        )
        h[Yindex(4,4,ℓₘᵢₙ)] += h̃₄₄
        h[Yindex(4,-4,ℓₘᵢₙ)] += -conj(h̃₄₄)
    end

    h
end
const mode_weights! = h!
