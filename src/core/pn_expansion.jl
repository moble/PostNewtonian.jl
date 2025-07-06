
@doc raw"""
    @pn_expansion pnsystem expansion

Mark `expansion` as a post-Newtonian expansion, to be truncated at the order given by the
`pnsystem`.

!!! warning

    The symbols `c`, `x`, `γₚₙ`, and `γ` are all considered to be PN-expansion terms inside
    the expression passed to this macro.  If present inside the expression, this macro will
    multiply each by the appropriate power of
    [`PNExpansionParameter(pnsystem)`](@ref PNExpansionParameter) to ensure that the
    expansion is truncated correctly.

Note that we generally interpret this PN order as a *relative* PN order.  For example, the
expression for [gravitational-wave energy flux](@ref 𝓕) looks like
```math
\mathcal{F} = \frac{32c^5}{5G} \nu^2 \left(\frac{v}{c}\right)^{10}
\left[1 - \left(\frac{ν}{12} + \frac{3}{4}\right) \left(\frac{v}{c}\right)^2 + \ldots\right]
```
Here, the part in square brackets is the expansion in terms of *relative* PN order, with the
1 representing the relative 0-pN term, the ``(v/c)^2`` term being the relative 1-pN term,
etc.  So in the code, we apply `@pn_expansion` to the part in square brackets:
```julia
32c^5 / 5G * ν^2 * (v / c)^10 * @pn_expansion(1 - (ν/12 + 3/4) * (v/c)^2 + ...)
```
Note that the entire expression is written inside a function that is modified by the
[`@pn_expression`](@ref) macro, which automatically inserts the `pnsystem` argument to this
`@pn_expansion` call, so you generally should not have to worry about it.

This macro gathers terms in `expansion` by the powers of ``1/c`` involved, zeroing out any
terms with powers of ``1/c`` higher than (twice) the `pnsystem`'s `PNOrder` parameter, and
combine the terms using the `PNExpansionReducer` specified in argument of the function that
includes this macro call.

The expansion and truncation are achieved by multiplying `c` by a
[`PNExpansionParameter(pnsystem)`](@ref PNExpansionParameter), which is just a
[`PNTerm`](@ref) with a coefficient of 1 and a `c⁻¹Exponent` of -1.  This redefinition
happens inside a `let` block created by this macro so that it doesn't interfere with any
`c`s on the outside.  For example, in the flux expression above, the `c`s in `32c^5 / 5G *
ν^2 * (v / c)^10` are outside the scope of the macro, so even if we set `c=2` above that
line the expansion will still be correctly truncated, and the value of `c` will be
preserved.

In the literature, some PN expansions do not *explicitly* use the `c` parameter, but instead
define the expansion in terms of the dimensionless parameters [`x`](@ref) or [`γₚₙ`](@ref),
which include powers of `c` in their definitions.  Thus, if these symbols are present inside
the expression passed to this macro, they will be multiplied by the appropriate power of
`PNExpansionParameter`.  The literature simply uses `γ`, rather than the more explicit `γₚₙ`
we use here, so `γ` will be treated the same as `γₚₙ`.

Euler's constant — sometimes called the Euler–Mascheroni constant, and distinct from Euler's
*number* ``e`` — is usually denoted by `γ` (including in Julia's own `Base.MathConstants`),
but is distinguished in the post-Newtonian literature by the subscript "E", and is defined
in this package as [`γₑ`](@ref).

!!! note

    Because `x` and `γₚₙ` are "of order" ``1/c^2``, their exponents are generally
    half-integers in any expansion.  Evaluating a term like `x^(7/2)` generically will be
    relatively inefficient, requiring a floating-point division, the logarithm of `x`, a
    floating-point multiplication, and exponentiation.  Compared to a single square-root and
    integer exponentiation, this is very slow.  Moreover, it does not work well with the
    [`PNTerm`](@ref) mechanism.  So another thing this macro does is search for terms like
    `x^(7/2)` and convert them to `(√x)^7`, and similarly for `γₚₙ` and `γ`.  This only
    works for explicit half-integers, rather than variables that may be half-integers.
"""
@public macro pn_expansion(pnsystem, expr)
    expansion_parameters = Expr[]
    if MacroTools.inexpr(expr, :c)
        push!(expansion_parameters, :(c = c * PNExpansionParameter($pnsystem)))
    end
    if MacroTools.inexpr(expr, :x)
        push!(expansion_parameters, :(x = x / PNExpansionParameter($pnsystem)^2))
        expr = MacroTools.postwalk(expr) do ex
            MacroTools.@capture(ex, x^(n_Int/2)) || return ex
            return :((√x)^$n)
        end
    end
    if MacroTools.inexpr(expr, :γₚₙ)
        push!(expansion_parameters, :(γₚₙ = γₚₙ / PNExpansionParameter($pnsystem)^2))
        expr = MacroTools.postwalk(expr) do ex
            MacroTools.@capture(ex, γₚₙ^(n_Int/2)) || return ex
            return :((√γₚₙ)^$n)
        end
    end
    if MacroTools.inexpr(expr, :γ)
        push!(expansion_parameters, :(γ = γ / PNExpansionParameter($pnsystem)^2))
        expr = MacroTools.postwalk(expr) do ex
            MacroTools.@capture(ex, γ^(n_Int/2)) || return ex
            return :((√γ)^$n)
        end
    end
    if isempty(expansion_parameters)
        error("No expansion parameters found in `@pn_expansion`.")
    end
    return esc(MacroTools.unblock(quote
        let $(expansion_parameters...)
            PNExpansionReducer($expr)
        end
    end))
end

@testitem "@pn_expansion" begin
    using PostNewtonian: @pn_expansion
    using MacroTools

    input = @macroexpand @pn_expansion pnsystem (
        1 - (ν/12 + 3/4) * (v/c)^2 + 4x^3 + 5x^(7/2)
    )
    output = MacroTools.unblock(
        quote
            let c = c * PNExpansionParameter(pnsystem),
                x = x / PNExpansionParameter(pnsystem)^2

                PNExpansionReducer(1 - (ν/12 + 3/4) * (v/c)^2 + 4x^3 + 5(√x)^7)
            end
        end,
    )
    @test MacroTools.striplines(input) == MacroTools.striplines(output)

    input = @macroexpand @pn_expansion pnsystem (1 - (ν/12 + 3/4)γ + 4γₚₙ^3 + 5γₚₙ^(7/2))
    output = MacroTools.unblock(
        quote
            let γₚₙ = γₚₙ / PNExpansionParameter(pnsystem)^2,
                γ = γ / PNExpansionParameter(pnsystem)^2

                PNExpansionReducer(1 - (ν/12 + 3/4)γ + 4γₚₙ^3 + 5(√γₚₙ)^7)
            end
        end,
    )
    @test MacroTools.striplines(input) == MacroTools.striplines(output)
end
