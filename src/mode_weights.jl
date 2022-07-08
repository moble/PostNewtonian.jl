# This function is stolen from SphericalFunctions.jl
@inline Yindex(ℓ, m, ℓₘᵢₙ=0) = ℓ*(ℓ + 1) - ℓₘᵢₙ^2 + m + 1


"""
    h!(h, u; ℓmin=0)
    h!(h, M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v; ℓmin=0)
    mode_weights!(h, u; ℓmin=0)
    mode_weights!(h, M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v; ℓmin=0)

Compute mode weights of gravitational waves emitted by `pn` system, modifying
`h` in place.

These modes are computed in the "co-orbital" frame, in which the larger object
lies on the positive ``x`` axis, the smaller lies on the negative ``x`` axis,
and the instantaneous angular velocity is in the positive ``z`` direction.

The modes are stored in `h` in order of increasing ``ℓ`` and increasing ``m``,
with ``m`` iterating fastest, all the way up to the highest available mode,
``(8,8)``.

Because gravitational waves have spin weight -2, the ``(ℓ,m)=(0,0)``,
``(1,-1)``, ``(1,0)``, and ``(1,1)`` modes are always 0.  By default, we assume
that these modes are nonetheless included in `h`.  If that is not the case, set
`ℓmin` to the smallest ``ℓ`` value that should be present in the output data —
`ℓmin=2` being the most reasonable alternative.

All non-spinning terms are taken from [Blanchet
(2014)](https://doi-org.proxy.library.cornell.edu/10.12942/lrr-2014-2).  The
1PN spin-orbit term is from Eq. (3.22d) of [Kidder
(1995)](https://link.aps.org/doi/10.1103/PhysRevD.52.821).  The 1.5PN
spin-orbit term is from Eq. (3.22f) of Kidder (1995) and Eq. (F15b) of [Will
and Wiseman (1996)](https://link.aps.org/doi/10.1103/PhysRevD.54.4813).  The
2PN spin-orbit term is from Eq. (4.13) of [Buonanno, Faye, Hinderer
(2013)](https://link.aps.org/doi/10.1103/PhysRevD.87.044009), while the 2PN
spin-spin term is from Eq. (4.15) of that reference.

"""
function h!(h, M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v; ℓmin=0)
    χ⃗₁ = QuatVec(χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ)
    χ⃗₂ = QuatVec(χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ)
    R = Quaternion(Rʷ, Rˣ, Rʸ, Rᶻ)
    fill!(h, 0)  # Set everything to 0 just to be safe
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), ℓ̂=ℓ̂(R), n̂=n̂(R), λ̂=λ̂(R), logv=log(v)
        T = typeof(logv)
        let π=T(π), γₑ=T(eulergamma), log2=T(log2), log3halves=T(log3halves), log5halves=T(log5halves), 𝒾=im
            Sℓ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            Σℓ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            Sn = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ n̂
            Σn = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ n̂
            Sλ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ λ̂
            Σλ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ λ̂
            S⃗₁ = χ⃗₁*M₁^2
            S⃗₂ = χ⃗₂*M₂^2
            S₁ℓ = S⃗₁ ⋅ ℓ̂
            S₁n = S⃗₁ ⋅ n̂
            S₁λ = S⃗₁ ⋅ λ̂
            S₂ℓ = S⃗₂ ⋅ ℓ̂
            S₂n = S⃗₂ ⋅ n̂
            S₂λ = S⃗₂ ⋅ λ̂

            c = 2ν * v^2 * √(16π/5)

            # ell=2
            h[Yindex(2,0,ℓmin)] = c * (-5/(14√T(6)))
            h[Yindex(2,1,ℓmin)] = c * (
                v^1 * (𝒾 * δ / 3)
                + v^3 * (𝒾 * δ * (-17 + 20ν) / 84)
                + v^4 * ((δ * (1 + 2𝒾*π + 4log2)) / 6)
                + v^5 * (𝒾 * δ * (-172 + ν * (-2036 + 237ν)) / 1512)
                + v^6 * (δ*(-34𝒾*π - 17*(1 + 4log2) + 2ν * (353 + 6𝒾*π + 12log2)) / 168)
            )
            h[Yindex(2,2,ℓmin)] = c * (
                1
                + v^2 * ((-107 + 55ν)/42)
                + v^3 * (2π)
                + v^4 * ((-2173 + ν * (-7483 + 2047ν)) / 1512)
                + v^5 * (-24𝒾*ν + ((-107 + 34ν)π) / 21)
                + v^6 * (
                    (27027409//646800) - 856γₑ/105 + (ν*(-834555 + ν*(-729396 + 114635ν))) / 99792
                    + 41ν * π^2 / 96 + (2π*(214𝒾 + 35π))/105 - 1712log2/105
                    - (856//105)*logv
                )
                + v^7 * ((-2𝒾 * ν * (-501655 + 24396ν) + 15*(-2173 + 2ν*(-2459 + 560ν))π) / 11340)
            )

            # ell=3
            h[Yindex(3,0,ℓmin)] = c * v^5 * (-2𝒾 * √T(6//7) * ν / 5)
            h[Yindex(3,1,ℓmin)] = c * (
                v^1 * (𝒾 * δ / 12√T(14))
                + v^3 * (-𝒾 * δ * (4 + ν) / 18√T(14))
                + v^4 * (δ * (7 + 5𝒾*π + 10log2) / 60√T(14))
                + v^5 * (-𝒾 * δ * (-607 + ν*(272 + 247ν)) / 2376√T(14))
                + v^6 * (
                    δ * (-5𝒾 * (16 + 7ν)π + 2*(-56 + ν - 5*(16 + 7ν)*log2)) / 360√T(14)
                )
                + v^7 * (
                    𝒾 * δ / 12√T(14) * (
                        (10753397//1513512) - 2log2*((212//105) + log2) - (26//21)*γₑ + (π^2/6)
                        - 2𝒾*π*((41//105) + log2) + (ν/8)*(-(1738843//19305) + (41//8)*π^2)
                        + (327059//30888)*ν^2 - (17525//15444)*ν^3
                        + logv * (-26//21)
                    )
                )
            )
            h[Yindex(3,2,ℓmin)] = c * (
                v^2 * (√T(5//7) * (1 - 3ν) / 3)
                + v^4 * ((-193 + (145 - 73ν)*5ν) / 54√T(35))
                + v^5 * ((-15𝒾 + 66𝒾*ν + 10π - 30π*ν) / 3√T(35))
                + v^6 * ((-1451 + (-17387 + 3*(33342 - 5341ν)ν)ν) / 2376√T(35))
            )
            h[Yindex(3,3,ℓmin)] = c * (
                v^1 * (-3𝒾 * √T(15//224) * δ)
                + v^3 * (-3𝒾 * √T(15//56) * δ * (-2 + ν))
                + v^4 * (√T(243//70) * δ * (-7 - 5𝒾*π + 10log3halves) / 4)
                + v^5 * (-𝒾 * √T(3//542080) * δ * (369 + (-3676 + 887ν)ν))
                + v^6 * (δ * (-3645𝒾 * (-8 + 3ν)π + (40824 - 96206ν + 7290*(-8 + 3ν)log3halves)) / 216√T(210))
                + v^7 * (
                    ((-3𝒾)/4) * √T(15//14) * δ * (
                        (19388147//280280) + (492//35)log3halves - 18*log3halves^2 - (78//7)γₑ + (3//2)π^2
                        + 6𝒾 * π * (-41//35 + 3log3halves)
                        + (-7055//3432 + 41//64 * π^2)ν - (318841//17160)ν^2 + (8237//2860)ν^3
                        + logv * (-(39//7) * 4log2)
                    )
                )
            )

            # ell=4
            h[Yindex(4,0,ℓmin)] = c * (-1 / 504√T(2))
            h[Yindex(4,1,ℓmin)] = c * (
                v^3 * (𝒾 * δ * (1 - 2ν) / 84√T(10))
                + v^5 * (-𝒾 * δ * (404 + (-1011 + 332ν)ν) / 11088√T(10))
                + v^6 * (δ * (64 - 1661ν - 30𝒾*(-1 + 2ν)π + 60*(1 - 2ν)log2) / 2520√T(10))
            )
            h[Yindex(4,2,ℓmin)] = c * (
                v^2 * (√T(5) * (1 - 3ν) / 63)
                + v^4 * ((-1311 + 5*(805 - 57ν)ν) / 4158√T(5))
                + v^5 * ((-21𝒾 + ν*(84𝒾 - 30π) + 10π) / 63√T(5))
                + v^6 * ((9342351 + 7ν*(-5460759 + 115ν*(34822 + 3363ν))) / 22702680√T(5))
            )
            h[Yindex(4,3,ℓmin)] = c * (
                v^3 * (9𝒾 * δ * (-1 + 2ν) / 4√T(70))
                + v^5 * (3𝒾 * δ * (468 + (-1267 + 524ν)ν) / 176√T(70))
                + v^6 * (δ * (-5184 + 16301ν + 2430𝒾*(-1 + 2ν)π + 4860*(1 - 2ν)log3halves) / 360√T(70))
            )
            h[Yindex(4,4,ℓmin)] = c * (
                v^2 * (8 * √T(5//7) * (-1 + 3ν) / 9)
                + v^4 * (4 * (1779 + 5ν*(-1273 + 525ν)) / 297√T(35))
                + v^5 * ((160*(-1 + 3ν)π + 𝒾*(336 - 1193ν + 320*(-1 + 3ν)log2)) / 9√T(35))
                + v^6 * ((-9618039 + 7ν*(9793071 + 5ν*(-3231338 + 678291ν))) / 405405√T(35))
            )

            # ell=5
            h[Yindex(5,1,ℓmin)] = c * (
                v^3 * (𝒾 * δ * (1 - 2ν) / 288√T(385))
                + v^5 * (-𝒾 * δ * (179 + 4*(-88 + ν)ν) / 11232√T(385))
                + v^6 * (δ * (181 - 70𝒾*(-1 + 2ν)π + 140log2 - 28ν*(313 + 10log2)) / 20160√T(385))
            )
            h[Yindex(5,2,ℓmin)] = c * (
                v^4 * ((2 + 10*(-1 + ν)ν) / 27√T(55))
                + v^6 * ((-3911 + 7ν*(3079 + 35ν*(-118 + 33ν))) / 12285√T(55))
            )
            h[Yindex(5,3,ℓmin)] = c * (
                v^3 * (9𝒾 * √T(3//110) * δ * (-1 + 2ν) / 32)
                + v^5 * (3𝒾 * √T(3//110) * δ * (207 + 8ν*(-58 + 11ν)) / 416)
                + v^6 * (δ * (-395847 + 1171828ν + 153090𝒾*(-1 + 2ν)π - 306180*(-1 + 2ν)*log3halves) / 60480√T(330))
            )
            h[Yindex(5,4,ℓmin)] = c * (
                v^4 * ((-32 - 160*(-1 + ν)ν) / 9√T(165))
                + v^6 * (16*(4451 - 7ν*(3619 + 5ν*(-1042 + 339ν))) / 4095√T(165))
            )
            h[Yindex(5,5,ℓmin)] = c * (
                v^3 * (-625𝒾 * δ * (-1 + 2ν) / 96√T(66))
                + v^5 * (-625𝒾 * δ * (263 + 16ν*(-43 + 16ν)) / 3744√T(66))
                + v^6 * (δ * (565625 - 1481676ν - 218750𝒾*(-1 + 2ν)π + 437500*(-1 + 2ν)log5halves) / 6720√T(66))
            )

            # ell=6
            h[Yindex(6,1,ℓmin)] = c * v^5 * (𝒾 * δ * (-1 + ν) * (-1 + 3 * ν) / 8316√T(26))
            h[Yindex(6,2,ℓmin)] = c * (
                v^4 * ((2 + (-1 + ν) * 10ν) / 297√T(65))
                + v^6 * ((-81 + (59 + (-64 + 7ν)ν) * 7ν) / 2079√T(65))
            )
            h[Yindex(6,3,ℓmin)] = c * v^5 * (-81𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 616√T(65))
            h[Yindex(6,4,ℓmin)] = c * (
                v^4 * (-128 * √T(2//39) * (1 + 5 * (-1 + ν)ν) / 495)
                + v^6 * (-64 * √T(2//39) * (-93 + 7 * (71 + (-88 + 19ν)ν)ν) / 3465)
            )
            h[Yindex(6,5,ℓmin)] = c * v^5 * (3125𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 504√T(429))
            h[Yindex(6,6,ℓmin)] = c * (
                v^4 * (54 * (1 + 5 * (-1 + ν)ν) / 5√T(143))
                + v^6 * (27 * (-113 + 7 * (91 + (-128 + 39ν)ν)ν) / 35√T(143))
            )

            # ell=7
            h[Yindex(7,1,ℓmin)] = c * v^5 * (𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 864864√T(2))
            h[Yindex(7,2,ℓmin)] = c * v^6 * ((1 - (-1 + ν)^2 * 7ν) / 3003√T(3))
            h[Yindex(7,3,ℓmin)] = c * v^5 * (-243𝒾 * √T(3//2) * δ * (-1 + ν) * (-1 + 3ν) / 160160)
            h[Yindex(7,4,ℓmin)] = c * v^6 * (128√T(2//33) * (-1 + (-1 + ν)^2 * 7ν) / 1365)
            h[Yindex(7,5,ℓmin)] = c * v^5 * (15625𝒾 * δ * (-1 + ν) * (-1 + 3ν) / 26208√T(66))
            h[Yindex(7,6,ℓmin)] = c * v^6 * (-81√T(3//143) * (-1 + (-1 + ν)^2 * 7ν) / 35)
            h[Yindex(7,7,ℓmin)] = c * v^5 * (-16807𝒾 * √T(7//858) * δ * (-1 + ν) * (-1 + 3ν) / 1440)

            # ell=8
            h[Yindex(8,2,ℓmin)] = c * v^6 * (-(-1 + (-1 + ν)^2 * 7ν) / 9009√T(85))
            h[Yindex(8,4,ℓmin)] = c * v^6 * ((128√T(2//187) * (-1 + (-1 + ν)^2 * 7ν)) / 4095)
            h[Yindex(8,6,ℓmin)] = c * v^6 * ((-243√T(3//17017) * (-1 + (-1 + ν)^2 * 7ν)) / 35)
            h[Yindex(8,8,ℓmin)] = c * v^6 * ((16384√T(2//85085) * (-1 + (-1 + ν)^2 * 7ν)) / 63)

            # Symmetric spin terms
            h[Yindex(2,0,ℓmin)] += c * v^4 * -((M₂*(S₁λ - S₁n) + M₁*(S₂λ - S₂n)) * (M₂*(S₁λ + S₁n) + M₁*(S₂λ + S₂n))) / (√T(6) * M^4 * ν^2)
            h[Yindex(2,1,ℓmin)] += c * v^2 * (𝒾 * Σℓ / 2M^2)
            h[Yindex(2,1,ℓmin)] += c * v^4 * (𝒾 * (-86*Sℓ*δ + Σℓ*(139ν - 79)) / 42M^2)
            h[Yindex(2,2,ℓmin)] += c * v^3 * (-(6Sℓ + 2Σℓ*δ) / 3M^2)
            h[Yindex(2,2,ℓmin)] += c * v^4 * (
                M₂^2 * (6S₁ℓ^2 + 5S₁λ^2 - 15𝒾*S₁λ*S₁n - 11S₁n^2)
                + M₁*M₂ * (12S₁ℓ*S₂ℓ + 10S₁λ*S₂λ - 15𝒾*S₁n*S₂λ - 15𝒾*S₁λ*S₂n - 22S₁n*S₂n)
                + M₁^2 * (6S₂ℓ^2 + 5S₂λ^2 - 15𝒾*S₂λ*S₂n - 11S₂n^2)
            ) / (6M^4 * ν^2)
            h[Yindex(3,1,ℓmin)] += c * v^4 * (√T(14)𝒾 * (Sℓ*δ - 5Σℓ*(3ν - 1)) / 336M^2)
            h[Yindex(3,2,ℓmin)] += c * v^3 * (2√T(35) * (Sℓ + Σℓ*δ) / 21M^2)
            h[Yindex(3,3,ℓmin)] += c * v^4 * (3√T(210)𝒾 * (7Sℓ*δ - 3Σℓ*(3ν - 1)) / 112M^2)
            h[Yindex(4,1,ℓmin)] += c * v^4 * (√T(10)𝒾 * (Sℓ*δ - 3Σℓ*ν + Σℓ) / 336M^2)
            h[Yindex(4,3,ℓmin)] += c * v^4 * (9√T(70)𝒾 * (-Sℓ*δ + 3Σℓ*ν - Σℓ) / 112M^2)

            # Symmetrize everything
            for ℓ in 2:8
                for m in 1:ℓ
                    h[Yindex(ℓ,-m,ℓmin)] = ifelse(isodd(ℓ), -1, 1) * conj(h[Yindex(ℓ,m,ℓmin)])
                end
            end

            # Anti-symmetric spin terms
            h̃_2_0 = c * (
                v^2 * (√T(6)𝒾 * Σn / 6M^2)
                + v^4 * (√T(6)𝒾 * (255Sn*δ - Σn*(506ν - 45)) / 126M^2)
            )
            h[Yindex(2,0,ℓmin)] += h̃_2_0
            h̃_2_1 = c * (
                v^3 * ((4𝒾*Sλ + 25*Sn + 4𝒾*Σλ*δ + 13*Σn*δ) / 6M^2)
                + v^4 * -3 * (M₂*S₁ℓ + M₁*S₂ℓ) * (M₂*S₁n + M₁*S₂n) / (2M^4 * ν^2)
            )
            h[Yindex(2,1,ℓmin)] += h̃_2_1
            h[Yindex(2,-1,ℓmin)] += -conj(h̃_2_1)
            h̃_2_2 = c * (
                v^2 * (-(Σλ + 𝒾*Σn) / 2M^2)
                + v^4 * ((19*Sλ*δ + 182𝒾*Sn*δ - 43*Σλ*ν + 5*Σλ - 280𝒾*Σn*ν + 98𝒾*Σn) / 84M^2)
            )
            h[Yindex(2,2,ℓmin)] += h̃_2_2
            h[Yindex(2,-2,ℓmin)] += -conj(h̃_2_2)
            h̃_3_0 = c * v^4 * (√T(42) * (-17Sλ*δ + Σλ*(35ν - 9)) / 168M^2)
            h[Yindex(3,0,ℓmin)] += h̃_3_0
            h̃_3_1 = c * v^3 * (√T(14) * (𝒾*Sλ + Sn + δ*(𝒾*Σλ + Σn)) / 21M^2)
            h[Yindex(3,1,ℓmin)] += h̃_3_1
            h[Yindex(3,-1,ℓmin)] += conj(h̃_3_1)
            h̃_3_2 = c * v^4 * (√T(35) * (-Σλ*(83ν - 17) + 4𝒾*Σn*(55ν - 13) + 25δ * (Sλ - 4𝒾*Sn)) / 168M^2)
            h[Yindex(3,2,ℓmin)] += h̃_3_2
            h[Yindex(3,-2,ℓmin)] += conj(h̃_3_2)
            h̃_3_3 = c * v^3 * (√T(210)𝒾 * (Sλ + 𝒾*Sn + δ*(Σλ + 𝒾*Σn)) / 21M^2)
            h[Yindex(3,3,ℓmin)] += h̃_3_3
            h[Yindex(3,-3,ℓmin)] += conj(h̃_3_3)
            h̃_4_0 = c * v^4 * (√T(2)𝒾 * (Sn*δ - 3Σn*ν + Σn) / 168M^2)
            h[Yindex(4,0,ℓmin)] += h̃_4_0
            h̃_4_2 = c * v^4 * (√T(5) * (-13Σλ*(3ν - 1) + 14𝒾*Σn*(3ν - 1) + δ*(13Sλ - 14𝒾*Sn)) / 168M^2)
            h[Yindex(4,2,ℓmin)] += h̃_4_2
            h[Yindex(4,-2,ℓmin)] += -conj(h̃_4_2)
            h̃_4_4 = c * v^4 * (9√T(35) * (-3Σλ*ν + Σλ - 𝒾*Σn*(3ν - 1) + δ*(Sλ + 𝒾*Sn)) / 56M^2)
            h[Yindex(4,4,ℓmin)] += h̃_4_4
            h[Yindex(4,-4,ℓmin)] += -conj(h̃_4_4)

        end
    end
    h
end
h!(h, u; ℓmin=0) = h!(h, u...; ℓmin=ℓmin)
const mode_weights! = h!
