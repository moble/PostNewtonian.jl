# Variable names refer to as-yet-unknown coefficients from v2 of Bini and Damour (2013b)
const a₆ᶜ¹ = 0
const a₆₅ᶜ¹ = 0
const a₇ˡⁿ¹ = 0
const a₇ᶜ¹ = 0

"""
    𝓔(pnsystem)
    binding_energy(pnsystem)

Compute the gravitational binding energy of a compact binary.

Note that this may not be as useful as its derivative, [`𝓔′`](@ref), which is used as part
of the right-hand side for orbital evolutions.

The nonspinning orbital binding energy is known through 4pN.  The expressions through 3.5pN
here come from Eq. (233) of [Blanchet (2014)](https://doi.org/10.12942/lrr-2014-2).

The 4pN term from Eq. (5.2d) of [Jaranowski and Schäfer](https://arxiv.org/abs/1303.3225v1)
is known exactly, now that the ``ν``-linear piece is given as Eq. (32) of [Bini and Damour
(2013a)](https://arxiv.org/abs/1305.4884v1).  The remaining terms are not known exactly, but
[Bini and Damour (2013b)](https://arxiv.org/abs/1312.2503v2) have derived some terms, though
there is incomplete information, which are noted as the constants in this code.

The spin-orbit terms in the energy are now complete to 4.0pN (the last term is zero).  These
terms come from Eq. (4.6) of [Bohé et al. (2012)](https://arxiv.org/abs/1212.5520v2).

The spin-squared terms (by which we mean both spin-spin and spin-orbit squared terms) in the
energy are known to 3pN order, and given in Eq. (3.33) of [Bohé et al.
(2015)](https://arxiv.org/abs/1501.01529).

The tidal-coupling terms come in to the binding energy at relative 5pN order, and are known
to 6pN order.  These terms come from Eq. (2.11) of [Vines et al.
(2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for
mass ratios, where ``χ₁ = m₁/m`` in their notation; in particular, ``χ`` is not a spin
parameter.  Also note that ``λ̂ = λ₂ v^{10}/(m₁+m₂)^5``, and we need to add the coupling
terms again with ``1 ↔ 2``.  Finally, note the normalization difference, where a different
overall factor is used, leading to a sign difference.
"""
@pn_expression function 𝓔(pnsystem)
    -1//2 * μ * v^2 * @pn_expansion(
        # Non-spinning terms; Eq. (233) of Blanchet (2014)
        1
        + v^2 * (-ν/12 - 3//4)
        + v^4 * (-ν^2/24 + 19ν/8 - 27//8)
        + v^6 * (-35ν^3/5184 - 155ν^2/96 + (34445//576 - 205π^2/96)ν - 675//64)

        # Eq. (5.2d) of Jaranowski and Schäfer
        + v^8 * (
            -3969//128 + 77ν^4/31104 + 301ν^3/1728 + (-498449//3456 + 3157π^2/576)ν^2
            + (-123671//5760 + 1792ln2/15 + 9037π^2/1536 + 896γₑ/15)ν
            + 2ln(v) * (448ν/15)
        )

        # Below are the incomplete terms from Eq. (74) of https://arxiv.org/abs/1312.2503v2
        + v^10 * (
            -45927//512 + ν^5/512 + 55ν^4/512 + (-1353π^2/256 + 69423//512)ν^3
            + (-21337π^2/1024 + 3a₆ᶜ¹ - 896ln2/5 - 448γₑ/5 + 893429//2880)ν^2
            + (-228916843//115200 - 9976γₑ/35 + 729ln3/7 - 23672ln2/35 + 126779π^2/512)ν
            + 2ln(v) * (-4988ν/35 - 656ν^2/5)
        )
        + v^11 * (10ν/3 * (13696π/525 + ν*a₆₅ᶜ¹))
        + v^12 * (
            -264627//1024
            + 2717ν^6/6718464
            + 5159ν^5/248832
            + (272855π^2/124416 - 20543435//373248)ν^4
            + (
                1232γₑ/27 + 6634243π^2/110592
                - 11a₆ᶜ¹/2 - 71700787//51840  + 2464ln2/27
            )ν^3
            + (
                113176680983//14515200 + 18491π^4/2304
                + 246004ln2/105 + 112772γₑ/105 + a₆ᶜ¹*11//2 + a₇ˡⁿ¹*2//3
                + a₇ᶜ¹*11//3 - 86017789π^2/110592 - 2673ln3/14
            )ν^2
            + (
                -389727504721//43545600 + 74888ln2/243 - 7128ln3/7
                - 30809603π^4/786432 - 3934568γₑ/8505 + 9118627045π^2/5308416
            )ν
            + 2ln(v) * (
                - 1967284ν/8505
                + 24464ν^3/135
                + (39754//105 + a₇ˡⁿ¹*11//3)ν^2
            )
        )

        # Spin-orbit; Eq. (4.6) of Bohé et al. (2012)
        + v^3 * (14sₗ/3 + 2δ * σₗ)
        + v^5 * ((11-61ν/9) * sₗ + (3 - 10ν/3)δ * σₗ)
        + v^7 * ((135//4 - 367ν/4 + 29ν^2/12) * sₗ + (27//4 - 39ν + 5ν^2/4)δ * σₗ)

        # Spin-squared; Eq. (3.33) of Bohé et al. (2015)
        + v^4 * (
            sₗ^2 * (-κ₊ - 2)
            + sₗ * σₗ * (-δ*κ₊ - 2δ + κ₋)
            + σₗ^2 * (δ*κ₋/2 - κ₊/2 + (κ₊ + 2)ν)
        )
        + v^6 * (
            sₗ^2 * (-5δ*κ₋/3 - 25*κ₊/6 + 50//9 + (5κ₊/6 + 5/3)ν)
            + sₗ * σₗ * (-5*δ*κ₊/2 + 25*δ/3 + 5κ₋/2 + (5δ*κ₊/6 + 5δ/3 + 35κ₋/6)ν)
            + σₗ^2 * (5δ*κ₋/4 - 5κ₊/4 + 5 + (5δ*κ₋/4 + 5κ₊/4 - 10)ν + (-5κ₊/6 - 5//3)ν^2)
        )

        # Spin-cubed; Eq. (6.17) of Marsat (2014)
        + v^7 * (
            sₗ^3 * (2κ₊ + 4λ₊ - 20)
            + sₗ^2 * σₗ * (2δ*κ₊ + 6δ*λ₊ - 32δ + 4κ₋ - 6λ₋)
            + sₗ * σₗ^2 * (5δ*κ₋ - 6δ*λ₋ - 5κ₊ + 6λ₊ - 12 + (-2κ₊ - 12λ₊ + 68)ν)
            + σₗ^3 * (-3δ*κ₊ + 2δ*λ₊ + 3κ₋ - 2λ₋ + (-2δ*λ₊ + 12δ - 6κ₋ + 6λ₋)ν)
        )

        # NS tidal coupling; Eq. (2.11) of Vines et al. (2011) with λ̂=v^10*Λ₂*(M₂/M)^5
        + v^10 * (
            - 9Λ₁ * ν * X₁^3
            - 9Λ₂ * ν * X₂^3
        )
        + v^12 * (
            - 11//2 * (3 + 2X₁ + 3X₁^2)Λ₁ * ν * X₁^3
            - 11//2 * (3 + 2X₂ + 3X₂^2)Λ₂ * ν * X₂^3
        )
    )
end
const binding_energy = 𝓔


# We derive the function 𝓔′ analytically from 𝓔.  Documentation goes below.
const 𝓔′ = let 𝓔=𝓔(symbolic_pnsystem), v=v(symbolic_pnsystem)
    ∂ᵥ = Differential(v)
    # Evaluate derivative symbolically
    𝓔′ = simplify(expand_derivatives(∂ᵥ(𝓔)), expand=true)#, simplify_fractions=false)
    # Turn it into (an Expr of) a function taking one argument: `pnsystem`
    𝓔′ = build_function(𝓔′, :pnsystem)
    # Remove `hold` (which we needed for Symbolics.jl to not collapse to Float64)
    𝓔′ = unhold(𝓔′)
    # "Flatten" the main sum, because Symbolics nests sums for some reason
    𝓔′ = apply_to_first_add!(𝓔′, flatten_add!)
    # Apply `@pn_expansion` to the main sum
    splitfunc = MacroTools.splitdef(𝓔′)
    splitfunc[:body] = apply_to_first_add!(
        splitfunc[:body],
        x->:(@pn_expansion(-1, $x))
    )
    𝓔′ = MacroTools.combinedef(splitfunc)
    # Finally, apply the "macro" to it and get a full function out
    eval(pn_expression(1, 𝓔′))::Function
end
const binding_energy_deriv=𝓔′

"""
    𝓔′(pnsystem)
    binding_energy_deriv(pnsystem)

Compute the derivative with respect to ``v`` of the binding energy of a compact binary.

This is computed symbolically from [`𝓔`](@ref); see that function for details.
"""
𝓔′
