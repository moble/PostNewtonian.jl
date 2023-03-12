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
    filter(n->getproperty(mod, n) isa T, names(mod, all=true))
end
