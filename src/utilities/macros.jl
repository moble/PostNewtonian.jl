"""
    type_converter(pnsystem, x)

Convert `x` to a type appropriate for the float type of `pnsystem`.

This is sort of an expansion of the `convert` function, but with nicer syntax for types from
this package, including the ability to do really weird things for `SymbolicPNSystem`.  It's
needed to ensure that the types of variables and constants are correct when we use them in
expressions, rather than just assuming everything is a `Float64`.
"""
function type_converter(pnsystem, x)
    return convert(eltype(pnsystem), x)
end
function type_converter(::FDPNSystem{FT}, x) where {FT}
    return x
end
function type_converter(::FDPNSystem{FT}, x::Integer) where {FT}
    return convert(FT, x)
end
function type_converter(::FDPNSystem{FT}, x::Rational) where {FT}
    return convert(FT, x)
end
function type_converter(::FDPNSystem{FT}, x::AbstractIrrational) where {FT}
    return convert(FT, x)
end

fundamental_variables = methodswith(PNSystem, FundamentalVariables)
fundamental_quaternionic_variables = [
    m for m ∈ methodswith(AbstractVector, FundamentalVariables) if m.name ∈ [:χ⃗₁, :χ⃗₂, :R]
]
derived_variables = methodswith(VecOrPNSystem, DerivedVariables)
pnvariables = map(v -> v.name, [fundamental_variables; derived_variables])

irrationals = unique(
    [
        find_symbols_of_type(Base.MathConstants, Irrational)
        find_symbols_of_type(MathConstants, Irrational)
    ],
)

# This should include all the unary functions that we want to use in any PN expression.
unary_funcs = [:√, :sqrt, :log, :ln, :sin, :cos]
# unary_funcs = Dict(
#     :√ => :(Base.sqrt),
#     :sqrt => :(Base.sqrt),
#     :log => :(Base.log),
#     :ln => :(Base.log),
# )

# for (expr, f) ∈ pairs(unary_funcs)
#     if expr ∈ (:sqrt, :ln)  # These are just aliases, so avoid redefinitions
#         continue
#     end
#     @eval begin
#         $f(::Type{T}, x) where T = $f(x)
#         $f(::Type{T}, x::Int) where T = $f(T(x))
#         $f(::Type{T}, x::Rational) where T = $f(T(x))
#         $f(::Type{T}, x::Irrational) where T = $f(T(x))
#     end
# end

function pn_expression(pnsystem::Symbol, body)
    # Look for variables in `body` that we need to treat specially, and write exprs to do
    # so.  These three are described as bullet points in the docstring of `@pn_expression`.
    pnvariables_exprs = [
        :($v = $v($pnsystem)) for v ∈ filter(v -> MacroTools.inexpr(body, v), pnvariables)
    ]
    irrationals_exprs = [
        :($v = type_converter($pnsystem, $v)) for
        v ∈ filter(v -> MacroTools.inexpr(body, v), irrationals)
    ]
    unary_funcs_exprs = [
        :($v = (x -> $v(type_converter($pnsystem, x)))) for
        v ∈ filter(v -> MacroTools.inexpr(body, v), unary_funcs)
    ]

    exprs = [
        pnvariables_exprs
        irrationals_exprs
        unary_funcs_exprs
    ]

    # Next, add pnsystem as the argument to each @pn_expansion call
    new_body = MacroTools.unblock(
        macroexpand(
            @__MODULE__,
            MacroTools.postwalk(body) do x
                if MacroTools.isexpr(x, :macrocall) &&
                    x.args[1] == Symbol("@pn_expansion") &&
                    !isa(x.args[end - 1], Symbol)
                    x′ = deepcopy(x)
                    insert!(x′.args, length(x′.args), pnsystem)
                    x′
                else
                    x
                end
            end;
            recursive=true,
        ),
    )

    # Finally, just wrap `new_body` in a `let` block, where we include exprs created above.
    # Also include the definitions `c=G=1` (to be overwritten inside any `@pn_expansion`).
    full_body = MacroTools.unblock(quote
        let c=G=1
            #@fastmath
            let $(exprs...)
                $(new_body)
            end
        end
    end)
    return full_body
end

function pn_expression(arg_index::Integer, func)
    splitfunc = MacroTools.splitdef(func)
    pnsystem = MacroTools.namify(splitfunc[:args][arg_index])
    splitfunc[:kwargs] = [
        splitfunc[:kwargs]
        :($(Expr(:kw, :(pn_expansion_reducer::Val{PNExpansionReducer}), :(Val(sum)))))
    ]
    splitfunc[:whereparams] = (splitfunc[:whereparams]..., :PNExpansionReducer)
    splitfunc[:body] = pn_expression(pnsystem, splitfunc[:body])
    return MacroTools.combinedef(splitfunc)
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
    return esc(pn_expression(1, func))
end

macro pn_expression(arg_index, func)
    return esc(pn_expression(arg_index, func))
end

"""
    @pn_expansion [pnsystem] expansion

Gather terms in `expansion` by the powers of ``1/c`` involved, zeroing out any terms with
powers of ``1/c`` higher than the `pnsystem`'s `PNOrder` parameter, and combine the terms
using the `PNExpansionReducer` specified in argument of the function that includes this
macro call.

This expansion is achieved by setting — inside a `let` block created by this macro —

Note that the `pnsystem` argument can be inserted automatically by [`@pn_expression`](@ref).
"""
macro pn_expansion(pnsystem, expr)
    return esc(pn_expansion(pnsystem, expr))
end

function pn_expansion(pnsystem, expr)
    quote
        let c = PNExpansionParameter($pnsystem)
            PNExpansionReducer($expr)
        end
    end
end
