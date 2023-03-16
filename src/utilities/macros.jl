non_pnvariables = [:eval, :include]
pnvariables = unique(filter(v->v∉non_pnvariables, [
    find_symbols_of_type(FundamentalVariables, Function);
    find_symbols_of_type(DerivedVariables, Function)
]))


for var ∈ pnvariables
    @eval begin
        function $var(v::PNSystem{T}) where {T<:Symbolics.Num}
            Symbolics.wrap(SymbolicUtils.Sym{Real}(Symbol($var)))
        end
    end
end

irrationals = unique([
    find_symbols_of_type(Base.MathConstants, Irrational);
    find_symbols_of_type(MathConstants, Irrational)
])

unary_funcs = [:√, :sqrt, :log, :ln, :sin, :cos]

function compute_pn_variables(arg_index, func)
    splitfunc = MacroTools.splitdef(func)
    pnsystem = splitfunc[:args][arg_index]
    body = splitfunc[:body]

    pnvariables_exprs = [
        :($v=PNVariables.$v($pnsystem))
        for v ∈ filter(v->MacroTools.inexpr(body, v), pnvariables)
    ]
    irrationals_exprs = [
        :($v=convert(eltype($pnsystem), $v))
        for v ∈ filter(v->MacroTools.inexpr(body, v), irrationals)
    ]
    unary_funcs_exprs = [
        :($v=($(esc(:x))->$v(convert(eltype($pnsystem), x))))
        for v ∈ filter(v->MacroTools.inexpr(body, v), unary_funcs)
    ]
    exprs = [
        pnvariables_exprs;
        irrationals_exprs;
        unary_funcs_exprs
    ]

    new_body = quote
        let $(exprs...)
            $body
        end
    end

    splitfunc[:body] = new_body
    MacroTools.combinedef(splitfunc)
end

"""
    @compute_pn_variables [arg_index=1] func

This macro takes the function `func`, looks for various symbols inside that function, and if
present defines them appropriately inside that function.

The first argument to this macro is `arg_index`, which just tells us which argument to the
function `func` is a `PNSystem`.  For example, the variables defined in
[`PostNewtonian.FundamentalVariables`](@ref Fundamental variables) all take a single
argument of `pnsystem`, which is used to compute the values for those variables.  The macro
just needs to know where to find `pnsystem`.

Once it has this information, there are three types of transformations it will make:

 1. For every [fundamental](@ref Fundamental variables) or [derived](@ref Derived variables)
    variable, the name of that variable used in the `func` will be replaced by its value
    when called with `pnsystem`.  For example, you can simply use [`M₁`](@ref) or
    [`μ`](@ref) in your code without have to call them as `M₁(pnsystem)` or `μ(pnsystem)`.
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
macro compute_pn_variables(func)
    compute_pn_variables(1, func)
end

macro compute_pn_variables(arg_index, func)
    compute_pn_variables(arg_index, func)
end
