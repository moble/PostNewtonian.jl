"""
    hold(x)

Delay evaluation of the argument in `Symbolics` expressions.

This is just a helper function that acts trivially — like the `identity` function — but also
gets registered with `Symbolics` to avoid evaluation of the argument.  For example, we can
preserve expressions like `π^2`, which Julia would normally convert directly to a `Float64`.

Note that you probably don't want to use this function directly; this will probably be done
for you by [`@pn_expression`](@ref) or similar.  If you *do* want to use this directly, you
probably want another layer of indirection to construct something like
`Symbolics.Num(SymbolicUtils.Term(hold, [x]))` so that you can use the result in a symbolic
expression.
"""
hold(x) = x
@register_symbolic hold(x)
Symbolics.derivative(::typeof(hold), args::NTuple{1,Any}, ::Val{1}) = 1

"""
    type_converter(pnsystem, x)

Convert `x` to a type appropriate for the float type of `pnsystem`.
"""
function type_converter(::PNSystem{T}, x) where {T<:Num}
    Symbolics.Num(SymbolicUtils.Term(hold, [x]))
end
function type_converter(::PNSystem{T}, x::Num) where {T<:Num}
    x
end
function type_converter(pnsystem, x)
    convert(eltype(pnsystem), x)
end

"""
    unhold(expr)

Remove occurrences of [`hold`](@ref) from an `Expr`.
"""
function unhold(expr)
    MacroTools.postwalk(expr) do x
        m = MacroTools.trymatch(:(f_(i_)), x)
        m === nothing || m[:f]!==hold ? x : Symbol(m[:i])
    end
end


fundamental_variables = methodswith(PNSystem, FundamentalVariables)
derived_variables = methodswith(PNSystem, DerivedVariables)
pnvariables = map(v->v.name, [fundamental_variables; derived_variables])

# Add symbolic capabilities to all derived variables (fundamental variables already work)
for method ∈ derived_variables
    name = method.name
    @eval begin
        function PostNewtonian.$name(v::PNSystem{T}) where {T<:Symbolics.Num}
            Symbolics.Num(SymbolicUtils.Sym{Real}(Symbol($name)))
        end
    end
end

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
    new_body = macroexpand(
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
    )

    # Finally, just wrap `new_body` in a `let` block, where we include exprs created above
    full_body = quote
        @fastmath let $(exprs...)
            $(new_body)
        end
    end

end

function pn_expression(arg_index::Integer, func)
    splitfunc = MacroTools.splitdef(func)
    pnsystem = MacroTools.namify(splitfunc[:args][arg_index])
    body = splitfunc[:body]
    splitfunc[:body] = pn_expression(pnsystem, body)
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

Once it has this information, there are four types of transformations it will make:

 1. For every [fundamental](@ref "Fundamental variables") or [derived](@ref "Derived
    variables") variable, the name of that variable used in the body of `func` will be
    replaced by its value when called with `pnsystem`.  For example, you can simply use the
    symbols `M₁` or `μ` in your code, rather than calling them as [`M₁(pnsystem)`](@ref M₁)
    or [`μ(pnsystem)`](@ref μ) every time they appear.
 2. Every `Irrational` defined in `Base.MathConstants` or `PostNewtonian.MathConstants` will
    be transformed to the `eltype` of `pnsystem`.  This lets you naturally use such
    constants in expressions like `2π/3` without automatically converting to `Float64`.
 3. Each of a short list of functions given by `unary_funcs` in `utilities/macros.jl` will
    first convert their arguments to the `eltype` of `pnsystem`.  In particular, you can use
    expressions like `√10` or `ln(2)` without the result being converted to a `Float64`.
 4. Insert the `pnsystem` argument as the first argument to each occurrence of
    `@pn_expansion` that needs it.

To be more explicit, the first three are achieved by defining the relevant quantities in a
`let` block placed around the body of `func`, so that the values may be used efficiently
without recomputation.
"""
macro pn_expression(func)
    esc(pn_expression(1, func))
end

macro pn_expression(arg_index, func)
    esc(pn_expression(arg_index, func))
end


"""
    @pn_expansion [pnsystem] expansion

Gather terms in `expansion` by the powers of `v` involved, the choose on the powers chosen
by the `pnsystem`'s `PNOrder` parameter, and evaluate efficiently in Horner form.

Note that the `pnsystem` argument can be inserted automatically by [`@pn_expression`](@ref).
For simplicity of presentation, we will assume that this is done in the examples below.

A "PN expansion" is a polynomial in ``v`` for which ``ln(v)`` factors may be present in
coefficients of ``v^k`` for ``k≥1``.  This function additionally requires that the input
should be a sum at its top level, meaning that this macro must be applied *after*
multiplication by any overall factor:
```julia
# This WILL NOT work:
@pn_expansion -M*ν*v^2/2 * (1 + v^2*(-ν/12-3//4))

# This WILL work:
-M*ν*v^2/2 * @pn_expansion(1 + v^2*(-ν/12-3//4))
```
Multiple terms may involve the same power of ``v``.  (I.e., the expansion does not have to
be simplified or
[`collect`](https://docs.sympy.org/latest/tutorials/intro-tutorial/simplification.html#collect)ed
in sympy parlance.)
"""
macro pn_expansion(pnsystem, expr)
    coefficients = var_collect(expr, :v)
    esc(:(evalpolysafe(v, $coefficients, $pnsystem)))
end

function evalpolysafe(v, coefficients, pnsystem)
    indices = 1:min(length(coefficients), order_index(pnsystem))
    evalpoly(v, coefficients[indices])
end


"""
    extract_var_factor(term, var)

Extract a factor of `var` from the product `term`.

This is a helper function for [`var_collect`](@ref).
"""
function extract_var_factor(term, var)
    if !MacroTools.isexpr(term, :call) || term.args[1] != :*
        if term == var
            return 1, 1
        end
        m = MacroTools.trymatch(:(v_^k_), term)
        if !isnothing(m) && m[:v] == var
            return m[:k], 1
        end
        return 0, term
    end
    term = deepcopy(term)
    k = 0
    indices = Int[]
    for (i,factor) ∈ enumerate(term.args)
        if i==1
            continue  # Skip the :*
        end
        if factor == var
            k += 1
            push!(indices, i)
            continue
        end
        m = MacroTools.trymatch(:(v_^k_), factor)
        if !isnothing(m) && m[:v] == var
            k += m[:k]
            push!(indices, i)
        end
    end
    splice!(term.args, indices)
    # TODO: If there's just 1 arg left over, return it alone
    if length(term.args) == 2
        k, term.args[2]
    else
        k, term
    end
end

"""
    var_collect(expr, var)

Collect coefficients in `expr` of various powers of the variable `var`.

The inputs should be an
[`Expr`](https://docs.julialang.org/en/v1/manual/metaprogramming/#Program-representation)
and a single `Symbol` to be found in that `Expr`.

The return value is an `Expr` representing a `Tuple` of values corresponding to the
coefficients of `var` to various powers.  For example,
```jl-doctest
julia> PostNewtonian.var_collect(:(1 + a*v + b*v^2 + c*v^4), :v)
:((1, a, b, 0, c))
```
(Note that there was *no* factor in `v^3`.)  This result is convenient for passing to
`evalpoly`, for example.
"""
function var_collect(expr, var)
    if !MacroTools.isexpr(expr, :call) || expr.args[1] != :+
        error("Input expression is not a sum at its highest level: $expr")
    end
    terms = Dict{Int,Any}()
    for term ∈ expr.args[2:end]
        k, term = extract_var_factor(term, var)
        if k ∈ keys(terms)
            terms[k] = :($(terms[k]) + $(term))
        else
            terms[k] = term
        end
    end
    term_exprs = [get(terms, k, 0) for k ∈ 0:maximum(keys(terms))]
    :(($(term_exprs...),))
end



# function symbolic_collect(expr::Num, var)
#     ex = Symbolics.unwrap(expr)
#     if Symbolics.operation(ex) != (+)
#         error("Input expression is not a sum")
#     end
#     terms = Dict{Int,Any}()
#     predicate(x) = isequal(x, var)
#     rulekLR = @acrule (~~a) * (~v::predicate)^(~k) * (~~b) => (~k, *(~~a..., ~~b...))
#     rulekL = @acrule (~~a) * (~v::predicate)^(~k) => (~k, *(~~a...))
#     rulekR = @acrule (~v::predicate)^(~k) * (~~b) => (~k, *(~~b...))
#     rulek = @acrule (~v::predicate)^(~k) => (~k, 1)
#     rule1LR = @acrule (~~a) * (~v::predicate) * (~~b) => (1, *(~~a..., ~~b...))
#     rule1L = @acrule (~~a) * (~v::predicate) => (1, *(~~a...))
#     rule1R = @acrule (~v::predicate) * (~~b) => (1, *(~~b...))
#     rule1 = @acrule (~v::predicate) => (1, 1)
#     rule0 = @acrule (~a::!predicate) => (0, a)
#     rules = SymbolicUtils.Chain([
#         rulekLR, rule1LR,
#         rulekL, rule1L,
#         rulekR, rule1R,
#         rulek, rule1,
#         rule0
#     ])
#     for x ∈ Symbolics.arguments(ex)
#         result = rules(x)
#         k, c = if result == x
#             0, x
#         else
#             result
#         end
#         if k ∈ keys(terms)
#             terms[k] = terms[k] + c
#         else
#             terms[k] = c
#         end
#     end
#     terms
# end
