"""
    𝓕(pnsystem)
    gw_energy_flux(pnsystem)

Compute the gravitational-wave energy flux to infinity

The nonspinning flux terms are complete to 4.5pN order, and are given in Eq. (6.11) of
[Blanchet et al. (2023)](https://arxiv.org/abs/2304.11186).

The spin-orbit terms in the flux are now known to 4.0pN.  These terms come from Eq. (4.9) of
[Marsat et al. (2013)](https://arxiv.org/abs/1307.6793v1)

The spin-squared terms (by which we mean both spin-spin and spin-orbit squared terms) in the
flux are known to 3pN order, and given most conveniently in Eq. (4.14) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).

Beyond 4.5pN, terms are only known in the extreme-mass-ratio limit.  These terms are given
in Appendix A of [Fujita (2012)](https://arxiv.org/abs/1211.5535v1).  He computed them up to
22pN.  That seems like overkill, so we'll just go up to 6pN.

For systems with matter, the tidal-heating terms come in at relative 5pN order, and are
known partially at 6pN order.  These terms come from Eq. (3.6) of [Vines et al.
(2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for
mass ratios, where ``χ₁ = m₁/m`` in their notation; in particular, ``χ`` is not a spin
parameter.  Also note that ``λ̂ = λ₂ v^{10}/(m₁+m₂)^5``, and we need to add the coupling
terms again with ``1 ↔ 2``.  Finally, note the normalization difference, where a different
overall factor is used, leading to a sign difference.
"""
@pn_expression function 𝓕(pnsystem)
    32ν^2/5 * v^10 * @pn_expansion(
        # Non-spinning terms; Eq. (314) of Blanchet (2014)
        1
        + v^2 * (-1247//336 - 35ν/12)
        + v^3 * (4π)
        + v^4 * (-44711//9072 + 9271ν/504 + 65ν^2/18)
        + v^5 * ((-8191//672 - 583*ν/24)π)
        + v^6 * (
            6643739519//69854400 + 16π^2/3 - 1712*(γₑ+2ln2+ln(v))/105
            + (-134543//7776 + 41π^2/48)ν - 94403ν^2/3024 - 775ν^3/324
        )
        + v^7 * ((-16285//504 + 214745ν/1728 + 193385ν^2/3024)π)
        + v^8 * (
            - 323105549467//3178375200 + 232597γₑ/4410 - 1369π^2/126
            + 39931ln2/294 - 47385ln3/1568 + 232597ln(v)/4410
            + (
                -1452202403629//1466942400 + 41478γₑ/245 - 267127π^2/4608
                + 479062ln2/2205 + 47385ln3/392 + 41478ln(v)/245
            )ν
            + (1607125//6804 - 3157π^2/384)ν^2 + 6875ν^3/504 + 5ν^4/6
        )
        + v^9 * (
            (
                265978667519//745113600 - 6848*(γₑ+2ln2+ln(v))/105
                + (2062241//22176 + 41π^2/12)ν - 133112905ν^2/290304 - 3719141ν^3/38016
            )π
        )

        # Spin-orbit terms; Eq. (4.9) of Marsat et al. (2013)
        + v^3 * (-4 * sₗ - 5δ/4 * σₗ)
        + v^5 * ((-9//2 + 272ν/9) * sₗ + (-13//16 + 43ν/4)δ * σₗ)
        + v^6 * ((-16π) * sₗ + (-31π/6)δ * σₗ)
        + v^7 * (
            (476645//6804 + 6172ν/189 - 2810ν^2/27) * sₗ
            + (9535//336 + 1849ν/126 - 1501ν^2/36)δ * σₗ
        )
        + v^8 * (
            (-3485//96 + 13879ν/72)π * sₗ
            + (-7163//672 + 130583ν/2016)π*δ * σₗ
        )

        # Spin-squared terms; Eq. (4.14) of Bohé et al. (2015)
        + v^4 * (
            sₗ^2 * (2κ₊ + 4)
            + sₗ * σₗ * (2δ*κ₊ + 4δ - 2κ₋)
            + σₗ^2 * (-δ*κ₋ + κ₊ + 1//16 + (-2κ₊ - 4)ν)
        )
        + v^6 * (
            sₗ^2 * (41δ*κ₋/16 - 271κ₊/112 - 5239//504 + (-43κ₊/4 - 43//2)ν)
            + sₗ * σₗ * (-279δ*κ₊/56 - 817δ/56 + 279κ₋/56 + (-43δ*κ₊/4 - 43δ/2 + κ₋/2)ν)
            + σₗ^2 * (
                279δ*κ₋/112 - 279*κ₊/112 - 25//8
                + (45δ*κ₋/16 + 243κ₊/112 + 344//21)ν
                + (43κ₊/4 + 43//2)ν^2
            )
        )

        # EMRI terms; Appendix A of Fujita (2012), with lower-order terms removed because
        # they have since been incorporated into non-EMRI terms above
        + v^10 * (
            - 2500861660823683//2831932303200 - 424223π^2/6804 - 83217611ln2/1122660
            + 916628467γₑ/7858620 + 47385ln3/196 + 916628467ln(v)/7858620
        )
        + v^11 * (
            - 142155π*ln3/784 + 8399309750401π/101708006400 + 177293γₑ*π/1176
            + 8521283π*ln2/17640 + 177293π*ln(v)/1176
        )
        + v^12 * (
            - 271272899815409ln2/157329572400
            - 54784π^2*ln2/315 - 246137536815857γₑ/157329572400 - 437114506833ln3/789268480
            - 256π^4/45 - 27392γₑ*π^2/315 - 27392ζ3/105 - 37744140625ln5/260941824
            + 1465472γₑ^2/11025 + 5861888γₑ*ln2/11025 + 5861888ln2^2/11025
            + 2067586193789233570693//602387400044430000 + 3803225263π^2/10478160
            + ln(v) * (
                - 246137536815857//157329572400 - 27392π^2/315
                + 2930944γₑ/11025 + 5861888ln2/11025
                + 1465472ln(v)/11025
            )
        )

        # # NS tidal heating; Eq. (3.6) of Vines et al. (2011)
        # + v^10 * (((12 - 18M / M₂)λ₂ + (12 - 18M / M₁)λ₁) / M^5)
        # + v^12 * (
        #     (
        #         (704 + 1803M₂/M - 4501*(M₂/M)^2 + 2170*(M₂/M)^3)λ₂ / (28M₂/M)
        #         + (704 + 1803M₁/M - 4501*(M₁/M)^2 + 2170*(M₁/M)^3)λ₁ / (28M₁/M)
        #     ) / M^5
        # )
    )
end
const gw_energy_flux = 𝓕
