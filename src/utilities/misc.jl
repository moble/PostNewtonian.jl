"""
    value(x)

Return `x` or the value wrapped by the `Dual` number `x`
"""
value(x) = hasproperty(x, :value) ? getproperty(x, :value) : x
