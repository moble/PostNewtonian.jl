"""
    ð(u)
    ð(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v)
    gw_energy_flux(u)
    gw_energy_flux(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v)

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

Beyond 3.5pN, terms other than the 4.0pN spin-orbit are only known in the
extreme-mass-ratio limit.  These terms are given in Appendix A of [Fujita
(2012)](https://arxiv.org/abs/1211.5535v1).  He computed them up to 22pN.
That seems like overkill, so we'll just go up to 6pN.

"""
function ð(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v)
    Ïââ = QuatVec(ÏââË£, ÏââÊ¸, Ïââá¶»)
    Ïââ = QuatVec(ÏââË£, ÏââÊ¸, Ïââá¶»)
    R = Quaternion(RÊ·, RË£, RÊ¸, Rá¶»)
    M = Mâ + Mâ
    let Î½=Î½(Mâ,Mâ), Î´=Î´(Mâ,Mâ), âÌ=âÌ(R), logv=log(v)
        let Ï=oftype(logv, Ï), Î³â=oftype(logv, eulergamma), Î¶3=oftype(logv, Î¶3)
            let log2=oftype(logv, log2), log3=oftype(logv, log3), log5=oftype(logv, log5)
                Sâ = S(Mâ,Mâ,Ïââ,Ïââ) â âÌ
                Î£â = Î£(Mâ,Mâ,Ïââ,Ïââ) â âÌ
                Ïââ = Ïâ(Mâ, Mâ, Ïââ, Ïââ)
                Ïââ = Ïâ(Mâ, Mâ, Ïââ, Ïââ)
                Ïââ = Ïââ â âÌ
                Ïââ = Ïââ â âÌ
                ÏâÂ² = Ïââ â Ïââ
                ÏâÂ² = Ïââ â Ïââ
                Ïââ = Ïââ â Ïââ

                32Î½^2/5 * v^10 * (
                    # Non-spinning terms; Eq. (314) of Blanchet (2014)
                    1
                    + v^2 * (-1247//336 - 35Î½/12)
                    + v^3 * (4Ï)
                    + v^4 * (-44711//9072 + 9271Î½/504 + 65Î½^2/18)
                    + v^5 * ((-8191//672 - 583*Î½/24)Ï)
                    + v^6 * (
                        6643739519//69854400 + 16Ï^2/3 - 1712*(Î³â+2log2+logv)/105
                        + (-134543//7776 + 41Ï^2/48)Î½ - 94403Î½^2/3024 - 775Î½^3/324
                    )
                    + v^7 * ((-16285//504 + 214745Î½/1728 + 193385Î½^2/3024)Ï)

                    # Spin-orbit terms; Eq. (4.9) of Marsat et al. (2013)
                    + v^3 * ((-4 * Sâ - 5Î´/4 * Î£â) / M^2)
                    + v^5 * (((-9//2 + 272Î½/9) * Sâ + (-13//16 + 43Î½/4)Î´ * Î£â) / M^2)
                    + v^6 * (((-16Ï) * Sâ + (-31Ï/6)Î´ * Î£â) / M^2)
                    + v^7 * (
                        (
                            (476645//6804 + 6172Î½/189 - 2810Î½^2/27) * Sâ
                            + (9535//336 + 1849Î½/126 - 1501Î½^2/36)Î´ * Î£â
                        ) / M^2
                    )
                    + v^8 * (
                        (
                            (-3485//96 + 13879Î½/72)Ï * Sâ
                            + (-7163//672 + 130583Î½/2016)Ï*Î´ * Î£â
                        ) / M^2
                    )

                    # Spin-squared terms; Eq. (C10) of Arun et al. (2008)
                    + v^4 * (
                        (287//96 + Î½/24) * (Ïââ)^2 - (89//96 + 7Î½/24) * (ÏâÂ² + 2Ïââ + ÏâÂ²) / 4
                        + (287//96 - 12Î½) * (Ïââ)^2 + (-89//96 + 4Î½) * (ÏâÂ² - 2Ïââ + ÏâÂ²) / 4
                        + 287Î´/48 * Ïââ * Ïââ
                        - 89Î´/48 * (ÏâÂ² - ÏâÂ²)/4
                    )

                    # EMRI terms; Appendix A of Fujita (2012)
                    + v^8 * (
                        -1369Ï^2/126 - 323105549467//3178375200 - 47385log3/1568
                        + 232597Î³â/4410 + 39931log2/294 + 232597logv/4410
                    )
                    + v^9 * (
                        -13696Ï*log2 / 105 - 6848Î³â*Ï/105 + 265978667519Ï/745113600 - 6848Ï*logv/105
                    )
                    + v^10 * (
                        - 2500861660823683//2831932303200 - 424223Ï^2/6804 - 83217611log2/1122660
                        + 916628467Î³â/7858620 + 47385log3/196 + 916628467logv/7858620
                    )
                    + v^11 * (
                        - 142155Ï*log3/784 + 8399309750401Ï/101708006400 + 177293Î³â*Ï/1176
                        + 8521283Ï*log2/17640 + 177293Ï*logv/1176
                    )
                    + v^12 * (
                        - 271272899815409log2/157329572400
                        - 54784Ï^2*log2/315 - 246137536815857Î³â/157329572400 - 437114506833log3/789268480 - 256Ï^4/45
                        - 27392Î³â*Ï^2/315 - 27392Î¶3/105 - 37744140625log5/260941824 + 1465472Î³â^2/11025
                        + 5861888Î³â*log2/11025 + 5861888log2^2/11025 + 2067586193789233570693//602387400044430000
                        + 3803225263Ï^2/10478160
                        + logv * (
                            - 246137536815857//157329572400 - 27392Ï^2/315
                            + 2930944Î³â/11025 + 5861888log2/11025
                            + 1465472logv/11025
                        )
                    )
                )
            end
        end
    end
end
ð(u) = ð(u...)
const gw_energy_flux = ð


"""
    ðNS(u, Î»â, Î»â)
    ðNS(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v, Î»â, Î»â)
    gw_energy_flux_NS(u, Î»â, Î»â)
    gw_energy_flux_NS(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v, Î»â, Î»â)

Compute tidal NS contribution to the gravitational-wave energy flux to infinity

For systems with matter, the tidal-coupling terms come in at relative 5pN
order, and are known partially at 6pN order.  These terms come from Eq. (3.6)
of [Vines et al. (2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).
Note their unusual convention for mass ratios, where ``Ïâ = mâ/m`` in their
notation; in particular, ``Ï`` is not a spin parameter.  Also note that ``Î»Ì =
Î»â v^{10}/(mâ+mâ)^5``, and we need to add the coupling terms again with ``1 â
2``.  Finally, note the normalization difference, where a different overall
factor is used, leading to a sign difference.

"""
function ðNS(Mâ, Mâ, ÏââË£, ÏââÊ¸, Ïââá¶», ÏââË£, ÏââÊ¸, Ïââá¶», RÊ·, RË£, RÊ¸, Rá¶», v, Î»â, Î»â)
    M = Mâ + Mâ

    32Î½^2/5 * v^10 * (
        # NS tides; Eq. (3.6) of Vines et al. (2011)
        v^10 * (((12 - 18M / Mâ)Î»â + (12 - 18M / Mâ)Î»â) / M^5)
        + v^12 * (
            (
                (704 + 1803Mâ/M - 4501*(Mâ/M)^2 + 2170*(Mâ/M)^3)Î»â / (28Mâ/M)
                + (704 + 1803Mâ/M - 4501*(Mâ/M)^2 + 2170*(Mâ/M)^3)Î»â / (28Mâ/M)
            ) / M^5
        )
    )
end
ðNS(u, Î»â, Î»â) = ðNS(u..., Î»â, Î»â)
const gw_energy_flux_NS = ðNS
