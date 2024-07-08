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
compute ``b_{i+1}`` from the condition that the term proportional to ``v^{i+1}`` in the
product of the series and its inverse must be zero.  This gives
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

function truncated_series_inverse!(b, a)
    @assert length(b) == length(a)
    n = length(a)
    @inbounds @fastmath if n > 0
        b[0+1] = inv(a[0+1])
    end
    @inbounds @fastmath for i ∈ 0:n-2
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
    @assert length(a) == length(b)
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
    @assert length(a) == length(b)
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


@doc raw"""
    truncated_series_ratio(a, b)

Evaluate the truncated ratio of the series `a` and `b`, evaluated at expansion value 1.
This is relevant when the expansion is not in the dynamic variable `v`, for example, but in
powers of ``1/c`` as in post-Newtonian expansions.  (That is, when the `v` dependence is
already include in the input coefficients.)
"""
function truncated_series_ratio(a::NTuple{N1,T1}, b::NTuple{N2,T2}) where {N1,N2,T1,T2}
    N = max(N1, N2)
    T = promote_type(T1, T2)
    if N2 == 0
        throw(DomainError("truncated_series_ratio(a,b): b must have at least one term"))
    elseif N1 == 0
        return zero(T)
    end
    b⁻¹ = MVector{N, T}(undef)

    @inbounds @fastmath begin
        b⁻¹[0+1] = inv(b[0+1])
        for i ∈ 0:N2-2
            b⁻¹[i+1+1] = -b⁻¹[0+1] * sum((b[j+1]*b⁻¹[i+1-j+1] for j ∈ 1:i+1), init=zero(T))
        end
        for i ∈ N2-1:N-2
            b⁻¹[i+1+1] = -b⁻¹[0+1] * sum((b[j+1]*b⁻¹[i+1-j+1] for j ∈ 1:N2-1), init=zero(T))
        end

        a╱b = zero(T)
        for i1 ∈ 1:N1
            a╱b += a[i1] * sum((b⁻¹[i2] for i2 ∈ 1:N-i1+1), init=zero(T))
        end
        a╱b
    end
end

@testitem "truncated_series_ratio(a,b)" begin
    using Random
    using DoubleFloats
    import PostNewtonian: truncated_series_inverse, truncated_series_ratio
    Random.seed!(123)
    for T ∈ [Float32, Float64, Double64]
        for N ∈ 1:20
            A = rand(T, N)
            A[1] = one(T) + rand(T) / 100
            a = Tuple(A)
            x = rand(T)
            ϵ = sum(a) * N * eps(T)
            for N1 ∈ 1:N
                expected = sum(truncated_series_inverse(a))
                unit = zeros(T, N1)
                unit[1] = one(T)
                @test truncated_series_ratio(Tuple(unit), a) ≈ expected rtol=ϵ
            end
            @test truncated_series_ratio(a, a) ≈ 1 rtol=ϵ
            @test truncated_series_ratio(a, x.*a) ≈ 1/x rtol=ϵ
            @test truncated_series_ratio(x.*a, a) ≈ x rtol=ϵ
        end
    end
end
