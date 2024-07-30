"""
    value(x)

Return `x` or the value wrapped by the `Dual` number `x`
"""
value(x) = hasproperty(x, :value) ? getproperty(x, :value) : x

"""
    find_symbols_of_type(mod, T)

Given a module `mod` (not just its name, but the actual imported module), find all objects
inside that module that are instances of the given type `T`.  The returned quantity is a
vector of `Symbol`s.
"""
function find_symbols_of_type(mod, T)
    return filter(n -> getproperty(mod, n) isa T, names(mod))
end

"""
    iscall(x, symbols)

Return `true` if the `Expr` `x` is a call to any element of `symbols`.
"""
iscall(x, symbols) = MacroTools.isexpr(x, :call) && x.args[1] ∈ symbols

"""
    isadd(x)

Return `true` if the `Expr` `x` is a call to `(+)` or `:+`.
"""
isadd(x) = iscall(x, ((+), :+))

"""
    ismul(x)

Return `true` if the `Expr` `x` is a call to `(*)` or `:*`.
"""
ismul(x) = iscall(x, ((*), :*))

"""
    flatten_binary!(expr, symbols)

Flatten nested binary operations — that is, apply associativity repeatedly.
"""
function flatten_binary!(expr, symbols)
    while iscall(expr, symbols) && any(x -> iscall(x, symbols), expr.args[2:end])
        args = expr.args[2:end]
        i₊ = findfirst(x -> iscall(x, symbols), args)
        args′ = [first(symbols); args[1:(i₊ - 1)]; args[i₊].args[2:end]; args[(i₊ + 1):end]]
        expr.args[:] = args′[1:length(expr.args)]
        append!(expr.args, args′[(1 + length(expr.args)):end])
    end
    return expr
end

flatten_add!(expr) = flatten_binary!(expr, ((+), :+))
flatten_mul!(expr) = flatten_binary!(expr, ((*), :*))

"""
    apply_to_first_add!(expr, func)

Apply `func` to the first sub-expression found in a "prewalk"-traversal of `expr` that
satisfies [`isadd`](@ref).  If `func` acts in place, so does this function.  In either case,
the expression should be returned.
"""
function apply_to_first_add!(expr, func)
    found_add = false
    MacroTools.prewalk(expr) do x
        if !found_add && isadd(x)
            found_add = true
            func(x)
        else
            x
        end
    end
end
