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

The spin-cubed terms are known to 3.5pN order, and come from Eq. (6.17) of [Marsat
(2014)](https://arxiv.org/abs/1411.4118).

The tidal-coupling terms come in to the binding energy at relative 5pN order, and are known
to 6pN order.  These terms come from Eq. (2.11) of [Vines et al.
(2011)](https://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for
mass ratios, where ``χ₁ = m₁/m`` in their notation; in particular, ``χ`` is not a spin
parameter.  Also note that ``λ̂ = λ₂ v^{10}/(m₁+m₂)^5``, and we need to add the coupling
terms again with ``1 ↔ 2``.  Finally, note the normalization difference, where a different
overall factor is used, leading to a sign difference.
"""
@pn_expression function 𝓔(pnsystem)
    return -μ * c^2 * (v / c)^2 / 2 * @pn_expansion(
        # Non-spinning terms; Eq. (233) of Blanchet (2014)
        1 +
            (v / c)^2 * (-ν / 12 - 3//4) +
            (v / c)^4 * (-ν^2 / 24 + 19ν / 8 - 27//8) +
            (v / c)^6 *
            (-35ν^3 / 5184 - 155ν^2 / 96 + (34445//576 - 205π^2 / 96)ν - 675//64)

            # Eq. (5.2d) of Jaranowski and Schäfer
            +
            (v / c)^8 * (
                -3969//128 +
                77ν^4 / 31104 +
                301ν^3 / 1728 +
                (-498449//3456 + 3157π^2 / 576)ν^2 +
                (-123671//5760 + 1792ln(2) / 15 + 9037π^2 / 1536 + 896γₑ / 15)ν +
                2ln(v) * (448ν / 15)
            )

            # Below are the incomplete terms from Eq. (74) of https://arxiv.org/abs/1312.2503v2
            +
            (v / c)^10 * (
                -45927//512 +
                ν^5 / 512 +
                55ν^4 / 512 +
                (-1353π^2 / 256 + 69423//512)ν^3 +
                (-21337π^2 / 1024 + 3a₆ᶜ¹ - 896ln(2) / 5 - 448γₑ / 5 + 893429//2880)ν^2 +
                (
                    -228916843//115200 - 9976γₑ / 35 + 729ln(3) / 7 - 23672ln(2) / 35 +
                    126779π^2 / 512
                )ν +
                2ln(v) * (-4988ν / 35 - 656ν^2 / 5)
            ) +
            (v / c)^11 * (10ν / 3 * (13696π / 525 + ν * a₆₅ᶜ¹)) +
            (v / c)^12 * (
                -264627//1024 +
                2717ν^6 / 6718464 +
                5159ν^5 / 248832 +
                (272855π^2 / 124416 - 20543435//373248)ν^4 +
                (
                    1232γₑ / 27 + 6634243π^2 / 110592 - 11a₆ᶜ¹ / 2 - 71700787//51840 +
                    2464ln(2) / 27
                )ν^3 +
                (
                    113176680983//14515200 +
                    18491π^4 / 2304 +
                    246004ln(2) / 105 +
                    112772γₑ / 105 +
                    a₆ᶜ¹ * 11//2 +
                    a₇ˡⁿ¹ * 2//3 +
                    a₇ᶜ¹ * 11//3 - 86017789π^2 / 110592 - 2673ln(3) / 14
                )ν^2 +
                (
                    -389727504721//43545600 + 74888ln(2) / 243 - 7128ln(3) / 7 -
                    30809603π^4 / 786432 - 3934568γₑ / 8505 + 9118627045π^2 / 5308416
                )ν +
                2ln(v) *
                (-1967284ν / 8505 + 24464ν^3 / 135 + (39754//105 + a₇ˡⁿ¹ * 11//3)ν^2)
            )

            # Spin-orbit; Eq. (4.6) of Bohé et al. (2012)
            +
            (v / c)^3 * (14sₗ / 3 + 2δ * σₗ) +
            (v / c)^5 * ((11 - 61ν / 9) * sₗ + (3 - 10ν / 3)δ * σₗ) +
            (v / c)^7 *
            ((135//4 - 367ν / 4 + 29ν^2 / 12) * sₗ + (27//4 - 39ν + 5ν^2 / 4)δ * σₗ)

            # Spin-squared; Eq. (3.33) of Bohé et al. (2015)
            +
            (v / c)^4 * (
                sₗ^2 * (-κ₊ - 2) +
                sₗ * σₗ * (-δ * κ₊ - 2δ + κ₋) +
                σₗ^2 * (δ * κ₋ / 2 - κ₊ / 2 + (κ₊ + 2)ν)
            ) +
            (v / c)^6 * (
                sₗ^2 * (-5δ * κ₋ / 3 - 25 * κ₊ / 6 + 50//9 + (5κ₊ / 6 + 5//3)ν) +
                sₗ *
                σₗ *
                (
                    -5 * δ * κ₊ / 2 +
                    25 * δ / 3 +
                    5κ₋ / 2 +
                    (5δ * κ₊ / 6 + 5δ / 3 + 35κ₋ / 6)ν
                ) +
                σₗ^2 * (
                    5δ * κ₋ / 4 - 5κ₊ / 4 +
                    5 +
                    (5δ * κ₋ / 4 + 5κ₊ / 4 - 10)ν +
                    (-5κ₊ / 6 - 5//3)ν^2
                )
            )

            # Spin-cubed; Eq. (6.17) of Marsat (2014)
            +
            (v / c)^7 * (
                sₗ^3 * (2κ₊ + 4λ₊ - 20) +
                sₗ^2 * σₗ * (2δ * κ₊ + 6δ * λ₊ - 32δ + 4κ₋ - 6λ₋) +
                sₗ * σₗ^2 * (5δ * κ₋ - 6δ * λ₋ - 5κ₊ + 6λ₊ - 12 + (-2κ₊ - 12λ₊ + 68)ν) +
                σₗ^3 * (-3δ * κ₊ + 2δ * λ₊ + 3κ₋ - 2λ₋ + (-2δ * λ₊ + 12δ - 6κ₋ + 6λ₋)ν)
            )

            # NS tidal coupling; Eq. (2.11) of Vines et al. (2011) with λ̂=v^10*Λ₂*(M₂/M)^5
            +
            (v / c)^10 * (-9Λ₁ * ν * X₁^3 - 9Λ₂ * ν * X₂^3) +
            (v / c)^12 * (
                -11//2 * (3 + 2X₁ + 3X₁^2)Λ₁ * ν * X₁^3 -
                11//2 * (3 + 2X₂ + 3X₂^2)Λ₂ * ν * X₂^3
            )
    )
end
const binding_energy = 𝓔

# NOTE: This is a helper function for the `@generated` function `𝓔′`; this function
# actually computes the code Expr to be generated.  This has been factored out to make it
# easier to generate different methods.  Specifically, we need to generate different code
# for `ForwardDiff.Dual`` numbers, which are only used in an extension to the core package.
# As such, the code relies on methods that cannot be defined yet, but generated functions
# "are only permitted to call functions that were defined *before* the definition of the
# generated function."  So we have to generate another method at a later time.  Therefore,
# we factor out this code to minimize duplication.
function 𝓔′code(
    ::Type{PN}, ::Type{Val{PNExpansionReducer}}, ::Type{ScalarType}, ::Type{FloatType}
) where {ST,PNOrder,PN<:PNSystem{ST,PNOrder},PNExpansionReducer,ScalarType,FloatType}
    # Create a `PNSystem` with `FastDifferentiation` (henceforth FD) variables, using the
    # same PNOrder as the input `pnsystem`.
    fdpnsystem = FDPNSystem(FloatType, PNOrder)

    # FD expects a single vector of variables, so we concatenate the state vector with the
    # two tidal-coupling parameters
    vars = FastDifferentiation.Node[fdpnsystem.state; Λ₁(fdpnsystem); Λ₂(fdpnsystem)]

    # Now we evaluate 𝓔 using the FD variables.  This will expand all derived variables in
    # terms of the fundamental variables, but FD will take care of evaluating those
    # efficiently via common subexpression elimination (CSE).
    𝓔formula = 𝓔(fdpnsystem; pn_expansion_reducer=Val(PNExpansionReducer))

    # Now we take the derivative of 𝓔 with respect to v.
    𝓔′ = SVector(FastDifferentiation.derivative(𝓔formula, v(fdpnsystem)))

    # Turn that into an Expr (FD insists on making it a function)
    in_place = true
    init_with_zeros = false
    𝓔′expr = FastDifferentiation.make_Expr(𝓔′, vars, in_place, init_with_zeros)

    # Now, we use `MacroTools` to get the body of the function.
    𝓔′body = MacroTools.unblock(MacroTools.splitdef(𝓔′expr)[:body])

    # At this point, the function is just a long series of statements inside an `@inbounds`
    # block, which we will want later, but first we need to extract them.
    MacroTools.@capture(𝓔′body, @inbounds begin
        𝓔′statements__
    end) || throw(
        ArgumentError(
            "\nNo @inbounds block found in 𝓔′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\n𝓔′($pnsystem)" *
            "\n```",
        ),
    )

    # The 𝓔′statements are mostly what we want, except that the last line is a return
    # statement.  We want that result, but we don't to return it yet; we want to wrap that
    # result, so we just get that returned quantity here.
    MacroTools.@capture(𝓔′statements[end], return 𝓔′return_) || throw(
        ArgumentError(
            "\nNo return statement found in 𝓔′ expression." *
            "\nSomething may have changed in FastDifferentiation." *
            "\nOpen an issue citing this Julia call:" *
            "\n```julia" *
            "\nusing PostNewtonian" *
            "\n𝓔′($pnsystem)" *
            "\n```",
        ),
    )
    𝓔′statements[end] = 𝓔′return

    if PNExpansionReducer ≡ identity
        # When `pn_expansion_reducer=Val(identity)` is passed, we return a PNExpansion
        NMax = Int(2PNOrder + 1)
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{$(length(𝓔′)),$(ScalarType)}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(𝓔′statements...)
            end
            return PNExpansion{$(length(𝓔′)),$(ScalarType),$NMax}(Tuple(result))
        end
    else
        # Otherwise, FD produces a 1-tuple, so we just extract the value from that.
        return quote
            input_variables = SVector(pnsystem)
            result = MVector{1,$(ScalarType)}(undef)
            result .= 0
            @fastmath @inbounds begin
                $(𝓔′statements...)
            end
            return result[1]
        end
    end
end

"""
    𝓔′(pnsystem)
    binding_energy_deriv(pnsystem)

Compute the derivative with respect to ``v`` of the binding energy of a compact binary.

This is computed automatically (via `FastDifferentiation`) from [`𝓔`](@ref); see that
function for details of the PN formulas.
"""
@generated function 𝓔′(
    pnsystem::PNSystem{ST,PNOrder}; pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum)
) where {ST,PNOrder,PNExpansionReducer}
    𝓔′code(pnsystem, pn_expansion_reducer, eltype(ST), eltype(ST))
end

const binding_energy_deriv = 𝓔′
