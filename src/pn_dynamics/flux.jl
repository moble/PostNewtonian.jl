"""
    𝓕(pn)
    gw_energy_flux(pn)

Compute the gravitational-wave energy flux to infinity

The nonspinning flux terms are complete to 3.5pN order.  These terms are given
by Eq. (314) of [Blanchet (2014)](https://doi.org/10.12942/lrr-2014-2).

The spin-squared terms (by which I mean both spin-spin and spin-orbit squared
terms) in the flux are known only at 2pN order (from [Kidder
(1995)](https://link.aps.org/doi/10.1103/PhysRevD.52.821) and [Will and Wiseman
(1996)](https://link.aps.org/doi/10.1103/PhysRevD.54.4813)).  They are most
conveniently given in Eq. (C10) of [Arun et
al. (2008)](https://arxiv.org/abs/0810.5336v3)

The spin-orbit terms in the flux are now known to 4.0pN.  These terms come from
Eq. (4.9) of [Marsat et al. (2013)](https://arxiv.org/abs/1307.6793v1)

"""
function 𝓕(pn)
    @unpack pn
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), ℓ̂=ℓ̂(R), π=oftype(v, π), γₑ=oftype(v, eulergamma)
        let log2=oftype(v, log2), logv=log(v)
            Sₗ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            Σₗ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            χ⃗ₐ = χₐ(M₁, M₂, χ⃗₁, χ⃗₂)
            χ⃗ₛ = χₛ(M₁, M₂, χ⃗₁, χ⃗₂)
            χₛₗ = χ⃗ₛ ⋅ ℓ̂
            χₐₗ = χ⃗ₐ ⋅ ℓ̂
            χ₁² = χ⃗₁ ⋅ χ⃗₁
            χ₂² = χ⃗₂ ⋅ χ⃗₂
            χ₁₂ = χ⃗₁ ⋅ χ⃗₂

            32ν^2/5 * v^10 * (
                # Non-spinning terms; Eq. (314) of Blanchet (2014)
                1
                + v^2 * (-1247//336 - 35ν/12)
                + v^3 * (4π)
                + v^4 * (-44711//9072 + 9271ν/504 + 65ν^2/18)
                + v^5 * ((-8191//672 - 583*ν/24)π)
                + v^6 * (
                    6643739519//69854400 + 16π^2/3 - 1712*(γₑ+2log2+logv)/105
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

            )
        end
    end
end
const gw_energy_flux = 𝓕


"""
    𝓕EMRI(pn)
    gw_energy_flux_EMRI(pn)

Compute the EMRI terms contributing to gravitational-wave energy flux to infinity

Beyond 3.5pN, the higher-order terms are only known in the extreme-mass-ratio
limit.  These terms are given in Appendix A of [Fujita
(2012)](https://arxiv.org/abs/1211.5535v1).  He computed them up to 22pN.  That
seems like overkill, so we'll just go up to 6pN.

"""
function 𝓕EMRI(pn)
    @unpack pn
    let ν=ν(M₁,M₂), ℓ̂=ℓ̂(R), π=oftype(v, π), γₑ=oftype(v, eulergamma)
        let log2=oftype(v, log2), log3=log(oftype(v, 3)), log5=log(oftype(v, 5)), ζ3=oftype(v, ζ3), logv=log(v)

            32ν^2/5 * v^10 * (
                # EMRI terms; Appendix A of Fujita (2012)
                + v^8 * (
                    -1369π^2/126 - 323105549467//3178375200 - 47385log3/1568
                    + 232597γₑ/4410 + 39931log2/294 + 232597logv/4410
                )
                + v^9 * (
                    -13696π*log2 / 105 - 6848γₑ*π/105 + 265978667519π/745113600 - 6848π*logv/105
                )
                + v^10 * (
                    - 2500861660823683//2831932303200 - 424223π^2/6804 - 83217611log2/1122660
                    + 916628467γₑ/7858620 + 47385log3/196 + 916628467logv/7858620
                )
                + v^11 * (
                    - 142155π*log3/784 + 8399309750401π/101708006400 + 177293γₑ*π/1176
                    + 8521283π*log2/17640 + 177293π*logv/1176
                )
                + v^12 * (
                    - 271272899815409log2/157329572400
                    - 54784π^2*log2/315 - 246137536815857γₑ/157329572400 - 437114506833log3/789268480 - 256π^4/45
                    - 27392γₑ*π^2/315 - 27392ζ3/105 - 37744140625log5/260941824 + 1465472γₑ^2/11025
                    + 5861888γₑ*log2/11025 + 5861888log2^2/11025 + 2067586193789233570693//602387400044430000
                    + 3803225263π^2/10478160
                    + logv * (
                        - 246137536815857//157329572400 - 27392π^2/315
                        + 2930944γₑ/11025 + 5861888log2/11025
                        + 1465472logv/11025
                    )
                )

            )
        end
    end
end
const gw_energy_flux_EMRI = 𝓕EMRI


"""
    𝓕NS(pn, λ₁, λ₂)
    gw_energy_flux_NS(pn, λ₁, λ₂)

Compute tidal NS contribution to the gravitational-wave energy flux to infinity

For systems with matter, the tidal-coupling terms come in at relative 5pN
order, and are known partially at 6pN order.  These terms come from Eq. (3.6)
of [Vines et al. (2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).
Note their unusual convention for mass ratios, where ``χ₁ = m₁/m`` in their
notation; in particular, ``χ`` is not a spin parameter.  Also note that ``λ̂ =
λ₂ v^{10}/(m₁+m₂)^5``, and we need to add the coupling terms again with ``1 ↔
2``.  Finally, note the normalization difference, where a different overall
factor is used, leading to a sign difference.

"""
function 𝓕NS(pn, λ₁, λ₂)
    @unpack pn
    M = M₁ + M₂

    32ν^2/5 * v^10 * (
        # NS tides; Eq. (3.6) of Vines et al. (2011)
        v^10 * (((12 - 18M / M₂)λ₂ + (12 - 18M / M₁)λ₁) / M^5)
        + v^12 * (
            (
                (704 + 1803M₂/M - 4501*(M₂/M)^2 + 2170*(M₂/M)^3)λ₂ / (28M₂/M)
                + (704 + 1803M₁/M - 4501*(M₁/M)^2 + 2170*(M₁/M)^3)λ₁ / (28M₁/M)
            ) / M^5
        )
    )
end
const gw_energy_flux_NS = 𝓕NS
