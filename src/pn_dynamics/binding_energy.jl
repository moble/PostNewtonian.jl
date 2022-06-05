# Equation numbers below refer to v1 of Bini and Damour (2013b)
const a_6__c1 = 0 # not yet known
const a_6__ln1 = -144//5  # coefficient of nu in Eq. (64)
const a_65__c1 = 0 # not yet known
const a_65__ln1 = 0 # not yet known
const a_7__ln1 = 0 # not yet known
const a_7__c1 = 0 # not yet known


"""
    𝓔(pn)
    binding_energy(pn)

Compute the binding energy of a compact binary

The nonspinning orbital binding energy is known through 4pN.  The expressions
through 3.5pN here come from Eq. (233) of [Blanchet
(2014)](https://doi.org/10.12942/lrr-2014-2).

The 4pN term from Eq. (5.2d) of [Jaranowski and
Schäfer](https://arxiv.org/abs/1303.3225v1) is known exactly, now that the
``ν``-linear piece is given as Eq. (32) of [Bini and Damour
(2013a)](https://arxiv.org/abs/1305.4884v1).  The remaining terms are not known
exactly, but [Bini and Damour (2013b)](https://arxiv.org/abs/1312.2503v1) have
derived some terms, though there is incomplete information, which are noted as
the constants in the following cell.  Note that, though the notation is
confusing, Bini and Damour claim they did not calculate the coefficient they
call ``a_6^{\\ln 1}``; but it seems to be given in their Eq. (64).

The spin-squared terms (by which I mean both spin-spin and spin-orbit squared
terms) in the energy are known only at 2pN order (from [Kidder
(1995)](https://link.aps.org/doi/10.1103/PhysRevD.52.821) and [Will and Wiseman
(1996)](https://link.aps.org/doi/10.1103/PhysRevD.54.4813)).  They are most
conveniently given in Eq. (C4) of [Arun et
al.](https://arxiv.org/abs/0810.5336v3)

The spin-orbit terms in the energy are now complete to 4.0pN (the last term is
zero).  These terms come from Eq. (4.6) of [Bohé et
al. (2012)](https://arxiv.org/abs/1212.5520v2):

"""
function 𝓔(pn)
    @unpack pn
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), ℓ̂=ℓ̂(R), logv=log(v)
        let γₑ=oftype(logv, eulergamma), π=oftype(logv, π), log2=oftype(logv, log2), log3=oftype(logv, log3)
            Sₗ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            Σₗ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            χ⃗ₐ = χₐ(M₁, M₂, χ⃗₁, χ⃗₂)
            χ⃗ₛ = χₛ(M₁, M₂, χ⃗₁, χ⃗₂)
            χₛₗ = χ⃗ₛ ⋅ ℓ̂
            χₐₗ = χ⃗ₐ ⋅ ℓ̂
            χ₁² = χ⃗₁ ⋅ χ⃗₁
            χ₂² = χ⃗₂ ⋅ χ⃗₂
            χ₁₂ = χ⃗₁ ⋅ χ⃗₂

            -M * ν * v^2 / 2 * (
                1
                + v^2 * (-3//4 - ν/12)
                + v^4 * (-27//8 + 19ν/8 - ν^2/24)
                + v^6 * (-675//64 + (34445//576 - 205π^2/96)ν - 155ν^2/96 - 35ν^3/5184)
                + v^8 * (
                    -3969//128 + (-123671//5760 + 9037π^2/1536 + 1792log2/15 + 896γₑ/15)ν
                    + (-498449//3456 + 3157π^2/576)ν^2 + 301ν^3/1728 + 77ν^4/31104
                    + logv * (896ν/15)
                )

                # Below are the incomplete terms
                + v^10 * (
                    -45927//512
                    + (-228916843//115200 - 9976γₑ/35 + 729log3/7 - 23672log2/35 + 126779π^2/512)ν
                    + (189745//576 + -21337π^2/1024 + 3a_6__c1 - 896log2/5 - 448γₑ/5 + 2a_6__ln1/3)ν^2
                    + (-1353π^2/256 + 69423//512)*ν^3
                    + 55ν^4/512
                    + ν^5/512
                    + logv * (-9976ν/35 + (-448//5 + 6a_6__ln1)ν^2)
                )
                + v^11 * (10ν/3 * (13696π/525 + ν*a_65__c1))
                + v^12 * (
                    -264627//1024
                    + (
                        -389727504721//43545600 + 74888log2/243 - 7128log3/7
                        - 3934568γₑ/8505 + 9118627045π^2/5308416 - 30809603π^4/786432
                    )ν 
                    + (
                        113594718743//14515200 + 18491π^4/2304
                        + 246004log2/105 + 112772γₑ/105 + 11a_6__c1/2 + a_6__ln1 + 2a_7__ln1/3
                        + 11a_7__c1/3 - 86017789π^2/110592 - 2673log3/14
                    )ν^2
                    + (
                        -75018547//51840 + 1232γₑ/27 + 6634243π^2/110592
                        - 11a_6__c1/2 + 2464log2/27 - 20a_6__ln1/9
                    )ν^3
                    + (272855π^2/124416 - 20543435//373248)ν^4
                    + 5159ν^5/248832
                    + 2717ν^6/6718464
                    + 2logv * (
                        11a_7__ln1/3
                        - 1967284ν/8505
                        + (56386//105 + 11a_6__ln1/2)*ν^2
                        + (616//27 - 11a_6__ln1/2)ν^3
                    )
                )

                # Spin-orbit
                + v^3 * ((14Sₗ/3 + 2δ * Σₗ) / M^2)
                + v^5 * (((11-61ν/9) * Sₗ + δ*(3 - 10ν/3) * Σₗ) / M^2)
                + v^7 * (((135//4 - 367ν/4 + 29ν^2/12) * Sₗ + δ*(27//4 - 39ν + 5ν^2/4) * Σₗ) / M^2)

                # Spin-squared
                + v^4 * (
                    (1 + δ - 2ν) * (χ₁² + χ₂²)/4 - 3*(χₐₗ^2 + χₛₗ^2)/2
                    - δ*(χ₂²/2 + 3χₐₗ*χₛₗ) + (χ₁₂ + 6χₐₗ^2)ν
                )

            )

        end
    end
end
const binding_energy = 𝓔


function 𝓔′(pn)
    @unpack pn
    M = M₁ + M₂
    let ν=ν(M₁,M₂), δ=δ(M₁,M₂), ℓ̂=ℓ̂(R), logv=log(v)
        let γₑ=oftype(logv, eulergamma), π=oftype(logv, π), log2=oftype(logv, log2), log3=oftype(logv, log3)
            Sₗ = S(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            Σₗ = Σ(M₁,M₂,χ⃗₁,χ⃗₂) ⋅ ℓ̂
            χ⃗ₐ = χₐ(M₁, M₂, χ⃗₁, χ⃗₂)
            χ⃗ₛ = χₛ(M₁, M₂, χ⃗₁, χ⃗₂)
            χₛₗ = χ⃗ₛ ⋅ ℓ̂
            χₐₗ = χ⃗ₐ ⋅ ℓ̂
            χ₁² = χ⃗₁ ⋅ χ⃗₁
            χ₂² = χ⃗₂ ⋅ χ⃗₂
            χ₁₂ = χ⃗₁ ⋅ χ⃗₂

            -M * ν * v / 2 * (
                2
                + 4v^2 * (-3//4 - ν/12)
                + 6v^4 * (-27//8 + 19ν/8 - ν^2/24)
                + 8v^6 * (-675//64 + (34445//576 - 205π^2/96)ν - 155ν^2/96 - 35ν^3/5184)
                + 10v^8 * (
                    -3969//128 + (-123671//5760 + 9037π^2/1536 + 1792log2/15 + 896γₑ/15)ν
                    + (-498449//3456 + 3157π^2/576)ν^2 + 301ν^3/1728 + 77ν^4/31104
                    + logv * (896ν/15)
                    + (896ν/150)
                )

                # Below are the incomplete terms
                + 12v^10 * (
                    -45927//512
                    + (-228916843//115200 - 9976γₑ/35 + 729log3/7 - 23672log2/35 + 126779π^2/512)ν
                    + (189745//576 + -21337π^2/1024 + 3a_6__c1 - 896log2/5 - 448γₑ/5 + 2a_6__ln1/3)ν^2
                    + (-1353π^2/256 + 69423//512)*ν^3
                    + 55ν^4/512
                    + ν^5/512
                    + logv * (-9976ν/35 + (-448//5 + 6a_6__ln1)ν^2)
                    + (-9976ν/35 + (-448//5 + 6a_6__ln1)ν^2)/12
                )
                + 13v^11 * (10ν/3 * (13696π/525 + ν*a_65__c1))
                + 14v^12 * (
                    -264627//1024
                    + (
                        -389727504721//43545600 + 74888log2/243 - 7128log3/7
                        - 3934568γₑ/8505 + 9118627045π^2/5308416 - 30809603π^4/786432
                    )ν 
                    + (
                        113594718743//14515200 + 18491π^4/2304
                        + 246004log2/105 + 112772γₑ/105 + 11a_6__c1/2 + a_6__ln1 + 2a_7__ln1/3
                        + 11a_7__c1/3 - 86017789π^2/110592 - 2673log3/14
                    )ν^2
                    + (
                        -75018547//51840 + 1232γₑ/27 + 6634243π^2/110592
                        - 11a_6__c1/2 + 2464log2/27 - 20a_6__ln1/9
                    )ν^3
                    + (272855π^2/124416 - 20543435//373248)ν^4
                    + 5159ν^5/248832
                    + 2717ν^6/6718464
                    + 2logv * (
                        11a_7__ln1/3
                        - 1967284ν/8505
                        + (56386//105 + 11a_6__ln1/2)*ν^2
                        + (616//27 - 11a_6__ln1/2)ν^3
                    )
                    + 2 * (
                        11a_7__ln1/3
                        - 1967284ν/8505
                        + (56386//105 + 11a_6__ln1/2)*ν^2
                        + (616//27 - 11a_6__ln1/2)ν^3
                    ) / 14
                )

                # Spin-orbit
                + 5v^3 * ((14Sₗ/3 + 2δ * Σₗ) / M^2)
                + 7v^5 * (((11-61ν/9) * Sₗ + δ*(3 - 10ν/3) * Σₗ) / M^2)
                + 9v^7 * (((135//4 - 367ν/4 + 29ν^2/12) * Sₗ + δ*(27//4 - 39ν + 5ν^2/4) * Σₗ) / M^2)

                # Spin-squared
                + 6v^4 * (
                    (1 + δ - 2ν) * (χ₁² + χ₂²)/4 - 3*(χₐₗ^2 + χₛₗ^2)/2
                    - δ*(χ₂²/2 + 3χₐₗ*χₛₗ) + (χ₁₂ + 6χₐₗ^2)ν
                )

            )

        end
    end
end
const binding_energy_deriv = 𝓔′


"""
    𝓔NS(pn)
    binding_energy_NS(pn)

Compute tidal NS contribution to the gravitational binding energy

The tidal-coupling terms come in to the energy at relative 5pN order, and are
known to 6pN order.  These terms come from Eq. (2.11) of [Vines et
al. (2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their
unusual convention for mass ratios, where ``χ₁ = m₁/m`` in their notation; in
particular, ``χ`` is not a spin parameter.  Also note that ``λ̂ = λ₂
v^{10}/(m₁+m₂)^5``, and we need to add the coupling terms again with ``1 ↔ 2``.
Finally, note the normalization difference, where a different overall factor is
used, leading to a sign difference.

"""
function 𝓔NS(pn)
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

            -M * ν * v^2 / 2 * (
                v^10 * (-9*((M₁/M₂)λ₂ + (M₂/M₁)λ₁) / M^5)
                + v^12 * (
                    (
                        -11//2*(M₁/M₂)*(3+2M₂/M+3*(M₂/M)^2)λ₂
                        - 11//2*(M₂/M₁)*(3+2M₁/M+3*(M₁/M)^2)λ₁
                    ) / M^5
                )
            )
        end
    end
end
const binding_energy_NS = 𝓔NS
