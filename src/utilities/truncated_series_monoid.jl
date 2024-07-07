@doc raw"""
    truncated_series_inverse(a)
    truncated_series_inverse!(b, a)

Given the coefficients `a` of a series, find the coefficients `b` of the multiplicative
inverse of that series, up to the order of the original series.

Note that this function returns the *coefficients* of the inverse, rather than its value.
This is relevant for use in [`truncated_series_product`](@ref) and
[`truncated_series_ratio`](@ref) — the latter of which just combines the former with this
function.

For example, suppose the input coefficients represent the series
```math
A \colonequals \sum_{i=0}^n a_{i-1} v^i.
```
(Remember that Julia's indexing is 1-based, so we subtract 1 to get the index of the
coefficient of ``v^i``.)  Then we return the coefficients `b` of the series
```math
B \colonequals \sum_{i=0}^n b_{i-1} v^i
```
such that
```math
A\, B = 1 + \mathcal{O}(v^{n+1}).
```

!!! note
    This function requires that `a[1]` be nonzero.  If you have a series that starts at
    a higher term — say, ``v^n`` for ``n>0`` — you should factor out the ``v^n``, and
    multiply the series resulting from this function by ``v^{-n}``.

## Explanation

The inverse coefficients can be computed fairly easily by induction.  Start by defining
```math
b_0 = 1/a_0.
```
Now, assuming that we've computed all coefficients up to and including ``b_{i}``, we can
compute
```math
b_{i+1} = -b_0\sum_{j=1}^{i} a_j b_{i-j}.
```
"""
function truncated_series_inverse(a::AbstractVector)
    b = similar(a)
    truncated_series_inverse!(b, a)
end

function truncated_series_inverse(a::NTuple{N, T}) where {N, T}
    b = MVector{N, T}(undef)
    truncated_series_inverse!(b, SVector{N, T}(a))
    Tuple(b)
end

@inbounds @fastmath function truncated_series_inverse!(b, a)
    n = length(a)
    if n > 0
        b[0+1] = inv(a[0+1])
    end
    for i ∈ 0:n-2
        b[i+1+1] = -b[0+1] * sum((a[j+1]*b[i+1-j+1] for j ∈ 1:i+1), init=zero(eltype(a)))
    end
    b
end

@doc raw"""
    truncated_series_product(a, b, v)

Evaluate the truncated product of the series `a` and `b`, which are expanded in powers of
`v`.

Note that this function returns the *value* of the summation, rather than its coefficients.

Here we define the series in terms of the coefficients `a` and `b` as
```math
A \colonequals \sum_{i=0}^n a_{i-1} v^i
\qquad
B \colonequals \sum_{i=0}^n b_{i-1} v^i,
```
and return the *value* of the product ``A\, B`` truncated at ``v^n``.

Internally, the sums are performed using `evalpoly`.

See also [`truncated_series_ratio`](@ref).
"""
function truncated_series_product(a, b, v)
    N = length(a)-1
    if N < 0
        return zero(v)
    end
    ex = b[N+1] * a[1]
    for n ∈ N-1:-1:0
        ex = muladd(v, ex, b[n+1] * evalpoly(v, @view a[1:(N-n)+1]))
    end
    ex
end

function truncated_series_product(a::NTuple{N,T}, b, v) where {N,T}
    if N < 1
        return zero(v)
    end
    ex = b[N] * a[1]
    for n ∈ N-2:-1:0
        ex = muladd(v, ex, b[n+1] * evalpoly(v, a[1:(N-n-1)+1]))
    end
    ex
end

@doc raw"""
    truncated_series_ratio(a, b, v)

Evaluate the truncated ratio of the series `a` and `b`, which are expanded in powers of `v`.

Note that this function returns the *value* of the summation, rather than its coefficients.

Here we define the series in terms of the coefficients `a` and `b` as
```math
A \colonequals \sum_{i=0}^n a_{i-1} v^i
\qquad
B \colonequals \sum_{i=0}^n b_{i-1} v^i,
```
and return the *value* of the ratio ``A / B`` truncated at ``v^n``.

This function simply combines [`truncated_series_product`](@ref) and
[`truncated_series_inverse`](@ref).
"""
function truncated_series_ratio(a, b, v)
    truncated_series_product(a, truncated_series_inverse(b), v)
end
