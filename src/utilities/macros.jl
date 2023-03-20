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

hold(x) = x
@register_symbolic hold(x)
Symbolics.derivative(::typeof(hold), args::NTuple{1,Any}, ::Val{1}) = 1
function type_converter(::PNSystem{T}, x) where {T<:Num}
    Symbolics.Num(SymbolicUtils.Term(hold, [x]))
end
function type_converter(::PNSystem{T}, x::Num) where {T<:Num}
    x
end
function type_converter(pnsystem, x)
    convert(eltype(pnsystem), x)
end
function unhold(expr)
    MacroTools.postwalk(expr) do x
        m = MacroTools.trymatch(:(f_(i_)), x)
        m === nothing || m[:f]!==hold ? x : Symbol(m[:i])
    end
end

function pn_expression(pnsystem::Symbol, body)
    irrationals_exprs = [
        :($v=type_converter($pnsystem, $v))
        for v ∈ filter(v->MacroTools.inexpr(body, v), irrationals)
    ]
    pnvariables_exprs = [
        :($v=$v($pnsystem))
        for v ∈ filter(v->MacroTools.inexpr(body, v), pnvariables)
    ]
    unary_funcs_exprs = [
        :($v=(x->$v(type_converter($pnsystem, x))))
        for v ∈ filter(v->MacroTools.inexpr(body, v), unary_funcs)
    ]
    exprs = [
        irrationals_exprs;
        pnvariables_exprs;
        unary_funcs_exprs
    ]

    new_body = quote
        #@fastmath let $(exprs...)
        let $(exprs...)
            $(body)
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
[`PostNewtonian.FundamentalVariables`](@ref Fundamental variables) all take a single
argument of `pnsystem`, which is used to compute the values for those variables; this macro
just needs to know where to find `pnsystem`.

Once it has this information, there are three types of transformations it will make:

 1. For every [fundamental](@ref Fundamental variables) or [derived](@ref Derived variables)
    variable, the name of that variable used in the body of `func` will be replaced by its
    value when called with `pnsystem`.  For example, you can simply use the symbols `M₁` or
    `μ` in your code, rather than calling them as [`M₁(pnsystem)`](@ref M₁) or
    [`μ(pnsystem)`](@ref μ) every time they appear.
 2. Every `Irrational` defined in `Base.MathConstants` or `PostNewtonian.MathConstants` will
    be transformed to the `eltype` of `pnsystem`.  This lets you naturally use such
    constants in expressions like `2π/3` without automatically converting to `Float64`.
 3. Each of a short list of functions given by `unary_funcs` in `utilities/macros.jl` will
    first convert their arguments to the `eltype` of `pnsystem`.  In particular, you can use
    expressions like `√10` or `ln(2)` without the result being converted to a `Float64`.

To be more precise, these are achieved by defining the relevant quantities in a `let` block
placed around the body of `func`, so that the values may be used efficiently without
recomputation.
"""
macro pn_expression(func)
    esc(pn_expression(1, func))
end

macro pn_expression(arg_index, func)
    esc(pn_expression(arg_index, func))
end
