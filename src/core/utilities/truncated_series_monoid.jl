@doc raw"""
    truncated_series_inverse(a)
    truncated_series_inverse!(b, a)

Given the coefficients `a` of a series, find the coefficients `b` of the *multiplicative*
inverse of that series, up to the order of the original series.  That is, if
```math
A \colonequals \sum_{i=0}^n a_{i} v^i,
```
then we return the coefficients `b` of the series
```math
B \colonequals \sum_{i=0}^n b_{i} v^i
```
such that
```math
A\, B = 1 + \mathcal{O}(v^{n+1}).
```

See [`lagrange_inversion`](@ref) for the *compositional* inverse (a.k.a. reversion), which
returns the coefficients of ``f^{-1}`` such that ``f^{-1}(f(v)) = v +
\mathcal{O}(v^{n+1})``.

Note that this function returns the *coefficients* of the inverse, rather than its value.
This is relevant for use in [`truncated_series_product`](@ref) and
[`truncated_series_ratio`](@ref) — the latter of which just combines the former with this
function.

!!! note

    This function requires that the constant term (`a[1]`) be nonzero.  If you have a series
    that starts at a higher term — say, ``v^n`` for ``n>0`` — you should factor out the
    ``v^n``, and multiply the series resulting from this function by ``v^{-n}``.

## Explanation

The inverse coefficients can be computed fairly easily by induction.  Start by defining
```math
b_0 = 1/a_0.
```
Now, assuming that we've computed all coefficients up to and including ``b_{i}``, we can
compute ``b_{i+1}`` from the condition that the term proportional to ``v^{i+1}`` in the
product of the series and its inverse must be zero.  That coefficient of that term is
clearly given by the sum of all pairs of coefficients ``a_j b_{i+1-j}`` for
``j=0,1,\ldots,i+1``:
```math
\sum_{j=0}^{i+1} a_j b_{i+1-j} = a_0 b_{i+1} + \sum_{j=1}^{i+1} a_j b_{i+1-j}.
```
Setting this last expression to zero, using the value of ``b_0`` above, and rearranging, we
have
```math
b_{i+1} = -b_0\sum_{j=1}^{i} a_j b_{i+1-j}.
```
"""
@public function truncated_series_inverse(a::AbstractVector)
    b = similar(a)
    return truncated_series_inverse!(b, a)
end

function truncated_series_inverse(a::NTuple{N,T}) where {N,T}
    b = MVector{N,T}(undef)
    truncated_series_inverse!(b, SVector{N,T}(a))
    return Tuple(b)
end

@public function truncated_series_inverse!(b, a)
    # We fake 0-based indexing by using `begin + i` for `i ∈ 0:n`.
    @assert length(b) == length(a)
    n = length(a) - 1
    @inbounds @fastmath if n ≥ 0
        b[begin + 0] = inv(a[begin + 0])
    end
    @inbounds @fastmath for i ∈ 0:(n - 1)
        b[begin + i + 1] =
            -b[begin + 0] * sum(
                (a[begin + j] * b[begin + i + 1 - j] for j ∈ 1:(i + 1));
                init=zero(eltype(a)),
            )
    end
    return b
end

@doc raw"""
    truncated_series_product(a, b, v)

Evaluate the truncated product of the series `a` and `b`, which are expanded in powers of
`v`.

Note that this function returns the *value* of the summation, rather than its coefficients.

Here we define the series in terms of the coefficients `a` and `b` as
```math
A \colonequals \sum_{i=0}^n a_{i} v^i
\qquad
B \colonequals \sum_{i=0}^n b_{i} v^i,
```
and return the *value* of the product ``A\, B`` truncated at ``v^n``.

Internally, the sums are performed using `evalpoly`.

See also [`truncated_series_ratio`](@ref).
"""
@public function truncated_series_product(a, b, v)
    @assert length(a) == length(b)
    N = length(a) - 1
    if N < 0
        return zero(v)
    end
    AB = b[begin + N] * a[begin + 0]
    for n ∈ (N - 1):-1:0
        AB = muladd(
            v, AB, b[begin + n] * evalpoly(v, @view a[(begin + 0):(begin + (N - n))])
        )
    end
    return AB
end

function truncated_series_product(a::NTuple{NT,T}, b, v) where {NT,T}
    @assert length(a) == length(b)
    N = length(a) - 1
    if N < 0
        return zero(v)
    end
    AB = b[begin + N] * a[begin + 0]
    for n ∈ (N - 1):-1:0
        AB = muladd(
            v, AB, b[begin + n] * evalpoly(v, @view a[(begin + 0):(begin + (N - n))])
        )
    end
    return AB
end

@doc raw"""
    truncated_series_ratio(a, b, v)

Evaluate the truncated ratio of the series `a` and `b`, which are expanded in powers of `v`.

Note that this function returns the *value* of the summation, rather than its coefficients.

Here we define the series in terms of the coefficients `a` and `b` as
```math
A \colonequals \sum_{i=0}^n a_{i} v^i
\qquad
B \colonequals \sum_{i=0}^n b_{i} v^i,
```
and return the *value* of the ratio ``A / B`` truncated at ``v^n``.

This function simply combines [`truncated_series_product`](@ref) and
[`truncated_series_inverse`](@ref).
"""
@public function truncated_series_ratio(a, b, v)
    return truncated_series_product(a, truncated_series_inverse(b), v)
end

@doc raw"""
    truncated_series_ratio(a, b)

Evaluate the truncated ratio of the series `a` and `b` at expansion value 1.

This is relevant when the expansion is not in the dynamic variable `v`, for example, but in
powers of ``1/c`` as in post-Newtonian expansions.  (That is, when the `v` dependence is
already included in the input coefficients.)

This is different from `truncated_series_ratio(a, b, 1)` in that `a` and `b` may have
different lengths, and it should be somewhat more efficient.
"""
function truncated_series_ratio(a::NTuple{N1,T1}, b::NTuple{N2,T2}) where {N1,N2,T1,T2}
    N = max(N1, N2)
    T = promote_type(T1, T2)
    if N2 == 0
        throw(DomainError("In truncated_series_ratio(a,b): b must have at least one term"))
    elseif N1 == 0
        return zero(T)
    end

    # The uppercase `N`s represent the number of terms in the tuples, while the lowercase
    # `n`s represent the highest index of the terms.
    n1, n2, n = N1 - 1, N2 - 1, N - 1

    # Next, we compute the same thing as `truncated_series_inverse` except that we truncate
    # this inverse series at `n`, instead of `n2`.  Specifically, the differences are that
    # (1) the range of iteration over `i` extends to `(n-1)` instead of `(n2-1)`, and (2)
    # the range of iteration over `j` extends to `min(i+1, n2)` instead of `(i+1)`.
    b⁻¹ = MVector{N,T}(undef)
    @inbounds @fastmath begin
        b⁻¹[begin + 0] = inv(b[begin + 0])
        for i ∈ 0:(n - 1)
            b⁻¹[begin + i + 1] =
                -b⁻¹[begin + 0] * sum(
                    (b[begin + j] * b⁻¹[begin + i + 1 - j] for j ∈ 1:min((i + 1), n2));
                    init=zero(T),
                )
        end

        # Now, we do the same thing as `truncated_series_product`, except that we account
        # for the fact that `a` may not be as long as `b⁻¹`, and we are assuming `v=1`.
        a╱b = zero(T)
        for i1 ∈ 0:n1
            a╱b += a[begin + i1] * sum((b⁻¹[begin + i2] for i2 ∈ 0:(n - i1)); init=zero(T))
        end
        a╱b
    end
end

@testitem "truncated_series_ratio(a,b)" begin
    using Random
    using DoubleFloats
    using StaticArrays: SVector
    import PostNewtonian: truncated_series_inverse, truncated_series_ratio
    Random.seed!(123)
    for T ∈ [Float32, Float64, Double64]
        for N ∈ 1:15
            A = rand(T, N)
            A[1] = one(T) + rand(T) / 100
            a = Tuple(A)
            x = rand(T)
            ϵ = eps(sum(a) * N)
            for N1 ∈ 1:N
                expected = sum(truncated_series_inverse(a))
                unit = zeros(T, N1)
                unit[1] = one(T)
                @test truncated_series_ratio(Tuple(unit), a) ≈ expected rtol = ϵ
            end
            @test truncated_series_ratio(a, a) ≈ 1 rtol = ϵ
            @test truncated_series_ratio(a, x .* a) ≈ 1 / x rtol = ϵ
            @test truncated_series_ratio(x .* a, a) ≈ x rtol = ϵ

            c = Tuple(rand(T, N))
            @test truncated_series_ratio(c, a) ≈
                truncated_series_ratio(SVector(c), SVector(a), 1) rtol = ϵ

            for Nzero ∈ 2:N
                cshort = c[begin:(Nzero - 1)]
                czeros = (cshort..., (zero(T) for i ∈ Nzero:N)...)
                #! format: off
                @test truncated_series_ratio(cshort, a) ≈
                    truncated_series_ratio(czeros, a) atol=2eps(T) rtol=2eps(T)

                ashort = a[begin:(Nzero - 1)]
                azeros = (ashort..., (zero(T) for i ∈ Nzero:N)...)
                @test truncated_series_ratio(c, ashort) ≈
                    truncated_series_ratio(c, azeros) atol=2eps(T) rtol=2eps(T)
                #! format: on
            end
        end
    end
end
