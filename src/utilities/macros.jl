

"""
    type_converter(pnsystem, x)

Convert `x` to a type appropriate for the float type of `pnsystem`.
"""
function type_converter(pnsystem, x)
    convert(eltype(pnsystem), x)
end
type_converter(pnsystem, x::FastDifferentiation.Node) = x
type_converter(::FDPNSystem{FT}, x::FastDifferentiation.Node) where FT = x

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

    # Finally, just wrap `new_body` in a `let` block, where we include exprs created above
    full_body = MacroTools.unblock(quote
        @fastmath let $(exprs...)
            $(new_body)
        end
    end)

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
    max_k, coefficients = var_collect(expr, :v)
    max_index = max_k + 1
    if offset != 0
        max_index_var = gensym("max_index")
        quote
            $max_index_var = min($max_index, order_index($pnsystem)-$offset)
            if $max_index_var < 1
                zero(v)
            else
                evalpoly(
                    v,
                    $coefficients[1:$max_index_var]
                )
            end
        end
    else
        :(evalpoly(
            v,
            $coefficients[1:min($max_index, order_index($pnsystem))]
        ))
    end
end


"""
    extract_var_factor(term, var)

Extract a factor of `var` from the product `term`.

This is a helper function for [`var_collect`](@ref).
"""
function extract_var_factor(term, var)
    if MacroTools.isexpr(term, :call) && term.args[1] ∈ ((/), :/)
        k₂, term₂ = extract_var_factor(term.args[2], var)
        k₃, term₃ = extract_var_factor(term.args[3], var)
        return k₂-k₃, Expr(:call, term.args[1], term₂, term₃)
        #return k₂-k₃, :($(term.args[1]), $term₂, $term₃)
    end
    if !MacroTools.isexpr(term, :call) || term.args[1] ∉ ((*), :*)
        if term == var
            return 1, 1
        end
        m = MacroTools.trymatch(:((^)(v_, k_)), term)
        m = !isnothing(m) ? m : MacroTools.trymatch(:($(^)(v_, k_)), term)
        if !isnothing(m) && m[:v] == var
            return m[:k], 1
        end
        return 0, term
    end
    term = flatten_mul!(deepcopy(term))
    k = 0
    indices = Int[]
    for (i,factor) ∈ enumerate(term.args)
        if i==1
            continue  # Skip the :*
        end
        if MacroTools.isexpr(factor, :call)
            k′, term′ = extract_var_factor(factor, var)
            # if term′ isa Expr
            #     term′ = Expr(:call, term′.args...)
            # end
            if k′ > 0
                k += k′
                term.args[i] = term′
            end
        else
            if factor == var
                k += 1
                push!(indices, i)
                continue
            end
            m = MacroTools.trymatch(:((^)(v_, k_)), factor)
            m = !isnothing(m) ? m : MacroTools.trymatch(:($(^)(v_, k_)), factor)
            if !isnothing(m) && m[:v] == var
                k += m[:k]
                push!(indices, i)
            end
        end
    end
    if !isempty(indices)
        splice!(term.args, indices)
    end
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

The return value is an `Int` representing the highest power of `v` in the expression, and an
`Expr` representing a `Tuple` of values corresponding to the coefficients of `var` to
various powers.  For example,
```jl-doctest
julia> PostNewtonian.var_collect(:(1 + a*v + b*v^2 + c*v^4), :v)
4, :((1, a, b, 0, c))
```
(Note that there was *no* factor in `v^3`.)  This result is convenient for passing to
`evalpoly`, for example.
"""
function var_collect(expr, var)
    if !MacroTools.isexpr(expr, :call)
        error("Input expression is not a call at its highest level: $expr")
    end
    terms = Dict{Int,Any}()
    if expr.args[1] ∉ ((+), :+, (-), :-)
        k, term = extract_var_factor(expr, var)
        terms[k] = term
    else
        for (i,term) ∈ enumerate(expr.args[2:end])
            k, term = extract_var_factor(term, var)
            if expr.args[1] ∈ ((-), :-) && i==2
                if k ∈ keys(terms)
                    terms[k] = :($(terms[k]) - $(term))
                else
                    terms[k] = :(-$(term))
                end
            else
                if k ∈ keys(terms)
                    terms[k] = :($(terms[k]) + $(term))
                else
                    terms[k] = term
                end
            end
        end
    end
    max_k = maximum(keys(terms))
    term_exprs = [get(terms, k, 0) for k ∈ 0:max_k]
    max_k, :(($(term_exprs...),))
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