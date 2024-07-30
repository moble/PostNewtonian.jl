"""

This file creates a very explicit binding energy function for a Post-Newtonian system.  It
largely duplicates the contents of the package's actual `binding_energy ≡ 𝓔` function, but
in a way that doesn't rely on macro fanciness or the PNTerm trickery.  This is useful for
testing the package's actual `binding_energy` function, as well as the validity of the
macro and PNTerm stuff.

Once the interface is finished, this could be done more elegantly with `@testsetup` from the
`TestItems` package.

"""

using PostNewtonian: a₆ᶜ¹, a₆₅ᶜ¹, a₇ˡⁿ¹, a₇ᶜ¹
function be(pnsystem, deriv)
    let M₁ = PostNewtonian.M₁(pnsystem),
        M₂ = PostNewtonian.M₂(pnsystem),
        v = PostNewtonian.v(pnsystem),
        Λ₁ = PostNewtonian.Λ₁(pnsystem),
        Λ₂ = PostNewtonian.Λ₂(pnsystem),
        X₁ = PostNewtonian.X₁(pnsystem),
        X₂ = PostNewtonian.X₂(pnsystem),
        M = PostNewtonian.M(pnsystem),
        sₗ = PostNewtonian.sₗ(pnsystem),
        δ = PostNewtonian.δ(pnsystem),
        μ = PostNewtonian.μ(pnsystem),
        ν = PostNewtonian.ν(pnsystem),
        σₗ = PostNewtonian.σₗ(pnsystem),
        χ₁² = PostNewtonian.χ₁²(pnsystem),
        χ₁₂ = PostNewtonian.χ₁₂(pnsystem),
        χ₂² = PostNewtonian.χ₂²(pnsystem),
        χₐₗ = PostNewtonian.χₐₗ(pnsystem),
        χₛₗ = PostNewtonian.χₛₗ(pnsystem),
        κ₊ = PostNewtonian.κ₊(pnsystem),
        κ₋ = PostNewtonian.κ₋(pnsystem),
        λ₊ = PostNewtonian.λ₊(pnsystem),
        λ₋ = PostNewtonian.λ₋(pnsystem),
        π = PostNewtonian.type_converter(pnsystem, π),
        ln2 = PostNewtonian.type_converter(pnsystem, PostNewtonian.ln2),
        ln3 = PostNewtonian.type_converter(pnsystem, PostNewtonian.ln3),
        γₑ = PostNewtonian.type_converter(pnsystem, PostNewtonian.γₑ),
        ln = (x -> log(PostNewtonian.type_converter(pnsystem, x))),
        pn_order = PostNewtonian.pn_order(pnsystem)

        e = Dict()
        eˡⁿ = Dict()

        c = -1//2 * μ  # NOTE: Included v^2 factor inside sum for easier differentiation
        e[0] = 1
        e[2] = (-ν / 12 - 3//4)
        e[4] = (-ν^2 / 24 + 19ν / 8 - 27//8)
        e[6] = (-35ν^3 / 5184 - 155ν^2 / 96 + (34445//576 - 205π^2 / 96)ν - 675//64)
        e[8] = (
            -3969//128 +
            77ν^4 / 31104 +
            301ν^3 / 1728 +
            (-498449//3456 + 3157π^2 / 576)ν^2 +
            (-123671//5760 + 1792ln2 / 15 + 9037π^2 / 1536 + 896γₑ / 15)ν
        )
        eˡⁿ[8] = (448ν / 15)

        # Below are the incomplete terms from Eq. (74) of https://arxiv.org/abs/1312.2503v2
        e[10] = (
            -45927//512 +
            ν^5 / 512 +
            55ν^4 / 512 +
            (-1353π^2 / 256 + 69423//512)ν^3 +
            (-21337π^2 / 1024 + 3a₆ᶜ¹ - 896ln2 / 5 - 448γₑ / 5 + 893429//2880)ν^2 +
            (
                -228916843//115200 - 9976γₑ / 35 + 729ln3 / 7 - 23672ln2 / 35 +
                126779π^2 / 512
            )ν
        )
        eˡⁿ[10] = (-4988ν / 35 - 656ν^2 / 5)
        e[11] = (10ν / 3 * (13696π / 525 + ν * a₆₅ᶜ¹))
        e[12] = (
            -264627//1024 +
            2717ν^6 / 6718464 +
            5159ν^5 / 248832 +
            (272855π^2 / 124416 - 20543435//373248)ν^4 +
            (
                1232γₑ / 27 + 6634243π^2 / 110592 - 11a₆ᶜ¹ / 2 - 71700787//51840 +
                2464ln2 / 27
            )ν^3 +
            (
                113176680983//14515200 +
                18491π^4 / 2304 +
                246004ln2 / 105 +
                112772γₑ / 105 +
                a₆ᶜ¹ * 11//2 +
                a₇ˡⁿ¹ * 2//3 +
                a₇ᶜ¹ * 11//3 - 86017789π^2 / 110592 - 2673ln3 / 14
            )ν^2 +
            (
                -389727504721//43545600 + 74888ln2 / 243 - 7128ln3 / 7 -
                30809603π^4 / 786432 - 3934568γₑ / 8505 + 9118627045π^2 / 5308416
            )ν
        )
        eˡⁿ[12] = (-1967284ν / 8505 + 24464ν^3 / 135 + (39754//105 + a₇ˡⁿ¹ * 11//3)ν^2)

        # Spin-orbit
        e[3] = (14sₗ / 3 + 2δ * σₗ)
        e[5] = ((11 - 61ν / 9) * sₗ + δ * (3 - 10ν / 3) * σₗ)
        e[7] = ((135//4 - 367ν / 4 + 29ν^2 / 12) * sₗ + δ * (27//4 - 39ν + 5ν^2 / 4) * σₗ)

        # Spin-squared
        e[4] += (
            sₗ^2 * (-κ₊ - 2) +
            sₗ * σₗ * (-δ * κ₊ - 2δ + κ₋) +
            σₗ^2 * (δ * κ₋ / 2 - κ₊ / 2 + (κ₊ + 2)ν)
        )
        e[6] += (
            sₗ^2 * (-5δ * κ₋ / 3 - 25 * κ₊ / 6 + 50//9 + (5κ₊ / 6 + 5//3)ν) +
            sₗ *
            σₗ *
            (-5 * δ * κ₊ / 2 + 25 * δ / 3 + 5κ₋ / 2 + (5δ * κ₊ / 6 + 5δ / 3 + 35κ₋ / 6)ν) +
            σₗ^2 * (
                5δ * κ₋ / 4 - 5κ₊ / 4 +
                5 +
                (5δ * κ₋ / 4 + 5κ₊ / 4 - 10)ν +
                (-5κ₊ / 6 - 5//3)ν^2
            )
        )

        # Spin-cubed
        e[7] += (
            sₗ^3 * (2κ₊ + 4λ₊ - 20) +
            sₗ^2 * σₗ * (2δ * κ₊ + 6δ * λ₊ - 32δ + 4κ₋ - 6λ₋) +
            sₗ * σₗ^2 * (5δ * κ₋ - 6δ * λ₋ - 5κ₊ + 6λ₊ - 12 + (-2κ₊ - 12λ₊ + 68)ν) +
            σₗ^3 * (-3δ * κ₊ + 2δ * λ₊ + 3κ₋ - 2λ₋ + (-2δ * λ₊ + 12δ - 6κ₋ + 6λ₋)ν)
        )

        # NS tidal coupling
        e[10] += -9 * (Λ₁ * X₁^3 + Λ₂ * X₂^3)ν
        e[12] += -11//2 * ((3 + 2X₁ + 3X₁^2)Λ₁ * X₁^3 + (3 + 2X₂ + 3X₂^2)Λ₂ * X₂^3)ν

        if deriv
            c *
            v *
            (
                sum(v^(k) * coeff * (k + 2) for (k, coeff) ∈ e if k ≤ 2pn_order; init=0) +
                sum(
                    v^(k) * coeff * 2 * ((k + 2) * log(v) + 1) for
                    (k, coeff) ∈ eˡⁿ if k ≤ 2pn_order;
                    init=0,
                )
            )
        else
            c *
            v^2 *
            (
                sum(v^(k) * coeff for (k, coeff) ∈ e if k ≤ 2pn_order; init=0) +
                sum(v^(k) * coeff * 2log(v) for (k, coeff) ∈ eˡⁿ if k ≤ 2pn_order; init=0)
            )
        end
    end
end
