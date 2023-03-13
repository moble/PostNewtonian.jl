"""
    𝓕(pnsystem)
    gw_energy_flux(pnsystem)

Compute the gravitational-wave energy flux to infinity

The nonspinning flux terms are complete to 3.5pN order.  These terms are given by Eq. (314)
of [Blanchet (2014)](https://doi.org/10.12942/lrr-2014-2).

The spin-squared terms (by which I mean both spin-spin and spin-orbit squared terms) in the
flux are known only at 2pN order (from [Kidder
(1995)](https://link.aps.org/doi/10.1103/PhysRevD.52.821) and [Will and Wiseman
(1996)](https://link.aps.org/doi/10.1103/PhysRevD.54.4813)).  They are most conveniently
given in Eq. (C10) of [Arun et al. (2008)](https://arxiv.org/abs/0810.5336v3)

The spin-orbit terms in the flux are now known to 4.0pN.  These terms come from Eq. (4.9) of
[Marsat et al. (2013)](https://arxiv.org/abs/1307.6793v1)

Beyond 3.5pN, terms other than the 4.0pN spin-orbit are only known in the extreme-mass-ratio
limit.  These terms are given in Appendix A of [Fujita
(2012)](https://arxiv.org/abs/1211.5535v1).  He computed them up to 22pN.  That seems like
overkill, so we'll just go up to 6pN.

For systems with matter, the tidal-coupling terms come in at relative 5pN order, and are
known partially at 6pN order.  These terms come from Eq. (3.6) of [Vines et al.
(2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for
mass ratios, where ``χ₁ = m₁/m`` in their notation; in particular, ``χ`` is not a spin
parameter.  Also note that ``λ̂ = λ₂ v^{10}/(m₁+m₂)^5``, and we need to add the coupling
terms again with ``1 ↔ 2``.  Finally, note the normalization difference, where a different
overall factor is used, leading to a sign difference.
"""
@compute_pn_variables function 𝓕(pnsystem)
    32ν^2/5 * v^10 * (
        # Non-spinning terms; Eq. (314) of Blanchet (2014)
        1
        + v^2 * (-1247//336 - 35ν/12)
        + v^3 * (4π)
        + v^4 * (-44711//9072 + 9271ν/504 + 65ν^2/18)
        + v^5 * ((-8191//672 - 583*ν/24)π)
        + v^6 * (
            6643739519//69854400 + 16π^2/3 - 1712*(γₑ+2ln2+lnv)/105
            + (-134543//7776 + 41π^2/48)ν - 94403ν^2/3024 - 775ν^3/324
        )
        + v^7 * ((-16285//504 + 214745ν/1728 + 193385ν^2/3024)π)

        # Spin-orbit terms; Eq. (4.9) of Marsat et al. (2013)
        + v^3 * ((-4 * Sₗ - 5δ/4 * Σₗ) / M^2)
        + v^5 * (((-9//2 + 272ν/9) * Sₗ + (-13//16 + 43ν/4)δ * Σₗ) / M^2)
        + v^6 * (((-16π) * Sₗ + (-31π/6)δ * Σₗ) / M^2)
        + v^7 * (
            (
                (476645//6804 + 6172ν/189 - 2810ν^2/27) * Sₗ
                + (9535//336 + 1849ν/126 - 1501ν^2/36)δ * Σₗ
            ) / M^2
        )
        + v^8 * (
            (
                (-3485//96 + 13879ν/72)π * Sₗ
                + (-7163//672 + 130583ν/2016)π*δ * Σₗ
            ) / M^2
        )

        # Spin-squared terms; Eq. (C10) of Arun et al. (2008)
        + v^4 * (
            (287//96 + ν/24) * (χₛₗ)^2 - (89//96 + 7ν/24) * (χ₁² + 2χ₁₂ + χ₂²) / 4
            + (287//96 - 12ν) * (χₐₗ)^2 + (-89//96 + 4ν) * (χ₁² - 2χ₁₂ + χ₂²) / 4
            + 287δ/48 * χₛₗ * χₐₗ
            - 89δ/48 * (χ₁² - χ₂²)/4
        )

        # EMRI terms; Appendix A of Fujita (2012)
        + v^8 * (
            -1369π^2/126 - 323105549467//3178375200 - 47385ln3/1568
            + 232597γₑ/4410 + 39931ln2/294 + 232597lnv/4410
        )
        + v^9 * (
            -13696π*ln2 / 105 - 6848γₑ*π/105 + 265978667519π/745113600 - 6848π*lnv/105
        )
        + v^10 * (
            - 2500861660823683//2831932303200 - 424223π^2/6804 - 83217611ln2/1122660
            + 916628467γₑ/7858620 + 47385ln3/196 + 916628467lnv/7858620
        )
        + v^11 * (
            - 142155π*ln3/784 + 8399309750401π/101708006400 + 177293γₑ*π/1176
            + 8521283π*ln2/17640 + 177293π*lnv/1176
        )
        + v^12 * (
            - 271272899815409ln2/157329572400
            - 54784π^2*ln2/315 - 246137536815857γₑ/157329572400 - 437114506833ln3/789268480
            - 256π^4/45 - 27392γₑ*π^2/315 - 27392ζ3/105 - 37744140625ln5/260941824
            + 1465472γₑ^2/11025 + 5861888γₑ*ln2/11025 + 5861888ln2^2/11025
            + 2067586193789233570693//602387400044430000 + 3803225263π^2/10478160
            + lnv * (
                - 246137536815857//157329572400 - 27392π^2/315
                + 2930944γₑ/11025 + 5861888ln2/11025
                + 1465472lnv/11025
            )
        )

        # NS tidal heating; Eq. (3.6) of Vines et al. (2011)
        + v^10 * (((12 - 18M / M₂)λ₂ + (12 - 18M / M₁)λ₁) / M^5)
        + v^12 * (
            (
                (704 + 1803M₂/M - 4501*(M₂/M)^2 + 2170*(M₂/M)^3)λ₂ / (28M₂/M)
                + (704 + 1803M₁/M - 4501*(M₁/M)^2 + 2170*(M₁/M)^3)λ₁ / (28M₁/M)
            ) / M^5
        )
    )
end
const gw_energy_flux = 𝓕
