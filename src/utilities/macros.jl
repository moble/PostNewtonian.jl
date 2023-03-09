irrationals = [
    find_symbols_of_type(Base.MathConstants, Irrational);
    find_symbols_of_type(MathConstants, Irrational)
]

pnvariables = filter(v->v!=:eval, [
    find_symbols_of_type(FundamentalVariables, Function);
    find_symbols_of_type(DerivedVariables, Function)
])

macro expand_variables(func)
    splitfunc = MacroTools.splitdef(func)
    arg1 = splitfunc[:args][1]
    body = splitfunc[:body]

    pnvariables_exprs = [
        :($v=PNVariables.$v($arg1))
        for v ∈ filter(v->MacroTools.inexpr(splitfunc[:body], v), pnvariables)
    ]
    irrationals_exprs = [
        :($v=oftype(eltype($arg1), $v))
        for v ∈ filter(v->MacroTools.inexpr(splitfunc[:body], v), irrationals)
    ]

    new_body = quote
        let $(pnvariables_exprs...)
            let $(irrationals_exprs...)
                $body
            end
        end
    end

    splitfunc[:body] = new_body
    MacroTools.combinedef(splitfunc)
end
