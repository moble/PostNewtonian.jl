irrationals = [
    find_symbols_of_type(Base.MathConstants, Irrational);
    find_symbols_of_type(MathConstants, Irrational)
]

pnvariables = filter(v->v!=:eval, [
    find_symbols_of_type(FundamentalVariables, Function);
    find_symbols_of_type(DerivedVariables, Function)
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
present defines them appropriately inside that function.  In particular, it defines PN
variables based on the value of an `PNSystem` argument to the function (located at
position `arg_index` in the argument list).  It also redefines `Irrational`s to have the
type relevant for that `PNSystem` object.
"""
macro compute_pn_variables(func)
    compute_pn_variables(1, func)
end

macro compute_pn_variables(arg_index, func)
    compute_pn_variables(arg_index, func)
end
