"""
    type_converter(pnsystem, x)

Convert `x` to a type appropriate for the float type of `pnsystem`.
"""
function type_converter(pnsystem, x)
    convert(eltype(pnsystem), x)
end
function type_converter(::FDPNSystem{FT}, x) where FT
    convert(FT, x)
end
function type_converter(::FDPNSystem, x::FastDifferentiation.Node)
    x
end

fundamental_variables = methodswith(PNSystem, FundamentalVariables)
fundamental_quaternionic_variables = [
    m for m ∈ methodswith(AbstractVector, FundamentalVariables)
    if m.name ∈ [:χ⃗₁, :χ⃗₂, :R]
]
derived_variables = methodswith(VecOrPNSystem, DerivedVariables)
pnvariables = map(v->v.name, [fundamental_variables; derived_variables])


irrationals = unique([
    find_symbols_of_type(Base.MathConstants, Irrational);
    find_symbols_of_type(MathConstants, Irrational)
])

unary_funcs = [:√, :sqrt, :log, :ln, :sin, :cos]


function pn_expression(pnsystem::Symbol, body)
    # Look for variables in `body` that we need to treat specially, and write exprs to do
    # so.  These three are described as bullet points in the docstring of `@pn_expression`.
    pnvariables_exprs = [
        :($v=$v($pnsystem))
        for v ∈ filter(v->MacroTools.inexpr(body, v), pnvariables)
    ]
    irrationals_exprs = [
        :($v=type_converter($pnsystem, $v))
        for v ∈ filter(v->MacroTools.inexpr(body, v), irrationals)
    ]
    unary_funcs_exprs = [
        :($v=(x->$v(type_converter($pnsystem, x))))
        for v ∈ filter(v->MacroTools.inexpr(body, v), unary_funcs)
    ]

    exprs = [
        pnvariables_exprs;
        irrationals_exprs;
        unary_funcs_exprs
    ]

    # Next, add pnsystem as the argument to each @pn_expansion call
    new_body = MacroTools.unblock(macroexpand(
        @__MODULE__,
        MacroTools.postwalk(body) do x
            if MacroTools.isexpr(x, :macrocall) &&
                x.args[1]==Symbol("@pn_expansion") &&
                !isa(x.args[end-1], Symbol)
                x′ = deepcopy(x)
                insert!(x′.args, length(x′.args), pnsystem)
                x′
            else
                x
            end
        end,
        recursive=true
    ))

    # Finally, just wrap `new_body` in a `let` block, where we include exprs created above.
    # Also include the definitions `c=G=1` (to be overwritten inside any `@pn_expansion`).
    full_body = MacroTools.unblock(quote
        c = one(eltype($pnsystem))
        G = one(eltype($pnsystem))
        @fastmath let $(exprs...)
            $(new_body)
        end
    end)

end

function pn_expression(arg_index::Integer, func)
    splitfunc = MacroTools.splitdef(func)
    pnsystem = MacroTools.namify(splitfunc[:args][arg_index])
    splitfunc[:kwargs] = [
        splitfunc[:kwargs];
        :($(Expr(:kw, :(pn_expansion_reducer::Val{PNExpansionReducer}), :(Val(sum)))))
    ]
    splitfunc[:whereparams] = (splitfunc[:whereparams]..., :PNExpansionReducer)
    splitfunc[:body] = pn_expression(pnsystem, splitfunc[:body])
    MacroTools.combinedef(splitfunc)
end

"""
    @pn_expression [arg_index=1] func

This macro takes the function `func`, looks for various symbols inside that function, and if
present defines them appropriately inside that function.

The first argument to this macro is `arg_index`, which just tells us which argument to the
function `func` is a `PNSystem`.  For example, the variables defined in
[`PostNewtonian.FundamentalVariables`](@ref "Fundamental variables") all take a single
argument of `pnsystem`, which is used to compute the values for those variables; this macro
just needs to know where to find `pnsystem`.

Once it has this information, there are five types of transformations it will make:

 1. Adds a keyword argument `pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum)` to the
    function signature.  This is used to determine how to reduce the PN expansion terms.
    The default is `Val(sum)`, which will just return a single number,  but `Val(identity)`
    can be used to return the expansion.  This should be used inside the function as
    `PNExpansionReducer`, and will be automatically used inside any `@pn_expansion`.
 2. For every [fundamental](@ref "Fundamental variables") or [derived](@ref "Derived
    variables") variable, the name of that variable used in the body of `func` will be
    replaced by its value when called with `pnsystem`.  For example, you can simply use the
    symbols `M₁` or `μ` in your code, rather than calling them as [`M₁(pnsystem)`](@ref M₁)
    or [`μ(pnsystem)`](@ref μ) every time they appear.
 3. Every `Irrational` defined in `Base.MathConstants` or `PostNewtonian.MathConstants` will
    be transformed to the `eltype` of `pnsystem`.  This lets you naturally use such
    constants in expressions like `2π/3` without automatically converting to `Float64`.
 4. Each of a short list of functions given by `unary_funcs` in `utilities/macros.jl` will
    first convert their arguments to the `eltype` of `pnsystem`.  In particular, you can use
    expressions like `√10` or `ln(2)` without the result being converted to a `Float64`.
 5. Insert the `pnsystem` argument as the first argument to each occurrence of
    `@pn_expansion` that needs it.

To be more explicit, the first three are achieved by defining the relevant quantities in a
`let` block placed around the body of `func`, so that the values may be used efficiently
without recomputation.

If you need to use one of the fundamental- or derived-variable functions as arguments of
values other than those encapsulated in `pnsystem`, you'll need to address them explicitly
with the module name — as in `PostNewtonian.v(;Ω, M)`.
"""
macro pn_expression(func)
    esc(pn_expression(1, func))
end

macro pn_expression(arg_index, func)
    esc(pn_expression(arg_index, func))
end


"""
    @pn_expansion [[offset] pnsystem] expansion

Gather terms in `expansion` by the powers of `v` involved, the choose on the powers chosen
by the `pnsystem`'s `PNOrder` parameter, and evaluate efficiently in Horner form.

Note that the `pnsystem` argument can be inserted automatically by [`@pn_expression`](@ref).
For simplicity of presentation, we will assume that this is done in the examples below.

A "PN expansion" is a polynomial in ``v`` for which ``ln(v)`` factors may be present in
coefficients of ``v^k`` for ``k≥1``.  The input may involve multiple terms with the same
power of ``v``.  (I.e., the expansion does not have to be simplified or
[`collect`](https://docs.sympy.org/latest/tutorials/intro-tutorial/simplification.html#collect)ed
in sympy parlance.)

However, note that the sum must appear together — as part of the same `:+` expression —
rather than, say, factored into two sums that multiply each other.

The `offset` argument represents the *index* of the relative PN order (that is, 2 times the
PN order) at which the PN expansion begins.  So if the first term in the expansion is really
a relative 2.5-pN term, we should use an `offset` of 5.  For example, the
[`tidal_heating`](@ref) expressions for `Ṁ₁` and `Ṁ₂` begin with a term proportional to
`v^15`, and are added to the flux [`𝓕`](@ref), which begin's with a term proportional to
`v^10`.  This means that `Ṁ₁` and `Ṁ₂` enter at 2.5-pN relative order, so we use an
`offset` of 5.
"""
macro pn_expansion(offset, pnsystem, expr)
    esc(pn_expansion(offset, pnsystem, expr))
end

macro pn_expansion(pnsystem, expr)
    esc(pn_expansion(0, pnsystem, expr))
end

function pn_expansion(offset, pnsystem, expr)
    quote
        let c = PNExpansionParameter($pnsystem)
            PNExpansionReducer($expr)
        end
    end
end
