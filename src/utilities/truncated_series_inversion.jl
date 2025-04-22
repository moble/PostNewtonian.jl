@doc raw"""
    x╱f_mod_xⁿ⁻¹(a)

Compute the truncated power series expansion of ``x/f`` mod ``x^{n-1}``, where ``f`` is the
power series expansion of a function ``f(x)``
```math
f(x) = a[1] x + a[2] x^2 + a[3] x^3 + \ldots + a[n-1] x^{n-1}.
```
Note, in particular that there is no constant term.  The result is a power series expansion
```math
h(x) = h[1] + h[2] x + h[3] x^2 + \ldots + h[n-1] x^{n-2},
```
which notably *does* have a constant term.

This function is essentially a helper function for the [`lagrange_inversion`](@ref)
function.

"""
function x╱f_mod_xⁿ⁻¹(a::NTuple{N,T}) where {N,T}
    b = zeros(MVector{N,typeof(inv(a[1]))})
    b[1] = inv(a[1])
    for i ∈ 2:N
        b[i] = -b[1] * sum((a[j] * b[i - j + 1] for j ∈ 2:i); init=zero(eltype(b)))
    end
    return Tuple(b)
end

@testitem "x╱f_mod_xⁿ⁻¹" begin
    using PostNewtonian: x╱f_mod_xⁿ⁻¹

    # We'll compare to exact values computed by SymPy for a few small values of `n`.
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{1}) = (1 / a[1],)
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{2}) = (1 / a[1], -a[2] / a[1]^2)
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{3}) = (
        1 / a[1], -a[2] / a[1]^2, (-a[1] * a[3] + a[2]^2) / a[1]^3
    )
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{4}) = (
        1 / a[1],
        -a[2] / a[1]^2,
        (-a[1] * a[3] + a[2]^2) / a[1]^3,
        (-a[1]^2 * a[4] + 2 * a[1] * a[2] * a[3] - a[2]^3) / a[1]^4,
    )
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{5}) = (
        1 / a[1],
        -a[2] / a[1]^2,
        (-a[1] * a[3] + a[2]^2) / a[1]^3,
        (-a[1]^2 * a[4] + 2 * a[1] * a[2] * a[3] - a[2]^3) / a[1]^4,
        (
            -a[1]^3 * a[5] + 2 * a[1]^2 * a[2] * a[4] + a[1]^2 * a[3]^2 -
            3 * a[1] * a[2]^2 * a[3] + a[2]^4
        ) / a[1]^5,
    )
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{6}) = (
        1 / a[1],
        -a[2] / a[1]^2,
        (-a[1] * a[3] + a[2]^2) / a[1]^3,
        (-a[1]^2 * a[4] + 2 * a[1] * a[2] * a[3] - a[2]^3) / a[1]^4,
        (
            -a[1]^3 * a[5] + 2 * a[1]^2 * a[2] * a[4] + a[1]^2 * a[3]^2 -
            3 * a[1] * a[2]^2 * a[3] + a[2]^4
        ) / a[1]^5,
        (
            -a[1]^4 * a[6] + 2 * a[1]^3 * a[2] * a[5] + 2 * a[1]^3 * a[3] * a[4] -
            3 * a[1]^2 * a[2]^2 * a[4] - 3 * a[1]^2 * a[2] * a[3]^2 +
            4 * a[1] * a[2]^3 * a[3] - a[2]^5
        ) / a[1]^6,
    )
    test_x╱f_mod_xⁿ⁻¹(a::NTuple{7}) = (
        1 / a[1],
        -a[2] / a[1]^2,
        (-a[1] * a[3] + a[2]^2) / a[1]^3,
        (-a[1]^2 * a[4] + 2 * a[1] * a[2] * a[3] - a[2]^3) / a[1]^4,
        (
            -a[1]^3 * a[5] + 2 * a[1]^2 * a[2] * a[4] + a[1]^2 * a[3]^2 -
            3 * a[1] * a[2]^2 * a[3] + a[2]^4
        ) / a[1]^5,
        (
            -a[1]^4 * a[6] + 2 * a[1]^3 * a[2] * a[5] + 2 * a[1]^3 * a[3] * a[4] -
            3 * a[1]^2 * a[2]^2 * a[4] - 3 * a[1]^2 * a[2] * a[3]^2 +
            4 * a[1] * a[2]^3 * a[3] - a[2]^5
        ) / a[1]^6,
        (
            -a[1]^5 * a[7] +
            2 * a[1]^4 * a[2] * a[6] +
            2 * a[1]^4 * a[3] * a[5] +
            a[1]^4 * a[4]^2 - 3 * a[1]^3 * a[2]^2 * a[5] - 6 * a[1]^3 * a[2] * a[3] * a[4] -
            a[1]^3 * a[3]^3 +
            4 * a[1]^2 * a[2]^3 * a[4] +
            6 * a[1]^2 * a[2]^2 * a[3]^2 - 5 * a[1] * a[2]^4 * a[3] + a[2]^6
        ) / a[1]^7,
    )

    for N ∈ 1:7
        a = (3, 7, 11, 17, 29, 37, 53)[1:N]
        @test all(test_x╱f_mod_xⁿ⁻¹(a) .≈ x╱f_mod_xⁿ⁻¹(a))
    end
end

@doc raw"""
    hⁱ✖h_mod_xⁿ⁻¹(h, i)

Compute the truncated power series expansion of ``h^i`` mod ``x^{n-1}``, where ``h`` is the
power series expansion of a function ``h(x)``
```math
h(x) = h[1] + h[2] x + h[3] x^2 + \ldots + h[n-1] x^{n-2}.
```
The result is a power series expansion
```math
h^i(x) = h^i[1] + h^i[2] x + h^i[3] x^2 + \ldots + h^i[n-1] x^{n-2}.
```

This function is essentially a helper function for the [`lagrange_inversion`](@ref)
function.

"""
function hⁱ✖h_mod_xⁿ⁻¹(hⁱ::NTuple{N,T}, h::NTuple{N,T}) where {N,T}
    hⁱ⁺¹ = zeros(MVector{N,T})
    for i ∈ 1:N
        hⁱ⁺¹[i] = sum((hⁱ[j] * h[i - j + 1] for j ∈ 1:i))
    end
    return Tuple(hⁱ⁺¹)
end

@testitem "hⁱ✖h_mod_xⁿ⁻¹" begin
    using PostNewtonian: hⁱ✖h_mod_xⁿ⁻¹

    # We'll compare to exact values computed by SymPy for a few small values of `n` and `i`
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{1}, ::Val{1}) = (h[0 + 1],)
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{1}, ::Val{2}) = (h[0 + 1]^2,)

    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{2}, ::Val{1}) = (h[0 + 1], h[1 + 1])
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{2}, ::Val{2}) = (h[0 + 1]^2, 2 * h[0 + 1] * h[1 + 1])
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{2}, ::Val{3}) = (h[0 + 1]^3, 3 * h[0 + 1]^2 * h[1 + 1])

    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{3}, ::Val{1}) = (h[0 + 1], h[1 + 1], h[2 + 1])
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{3}, ::Val{2}) = (
        h[0 + 1]^2, 2 * h[0 + 1] * h[1 + 1], 2 * h[0 + 1] * h[2 + 1] + h[1 + 1]^2
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{3}, ::Val{3}) = (
        h[0 + 1]^3,
        3 * h[0 + 1]^2 * h[1 + 1],
        3 * h[0 + 1]^2 * h[2 + 1] + 3 * h[0 + 1] * h[1 + 1]^2,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{3}, ::Val{4}) = (
        h[0 + 1]^4,
        4 * h[0 + 1]^3 * h[1 + 1],
        4 * h[0 + 1]^3 * h[2 + 1] + 6 * h[0 + 1]^2 * h[1 + 1]^2,
    )

    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{4}, ::Val{1}) = (h[0 + 1], h[1 + 1], h[2 + 1], h[3 + 1])
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{4}, ::Val{2}) = (
        h[0 + 1]^2,
        2 * h[0 + 1] * h[1 + 1],
        2 * h[0 + 1] * h[2 + 1] + h[1 + 1]^2,
        2 * h[0 + 1] * h[3 + 1] + 2 * h[1 + 1] * h[2 + 1],
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{4}, ::Val{3}) = (
        h[0 + 1]^3,
        3 * h[0 + 1]^2 * h[1 + 1],
        3 * h[0 + 1]^2 * h[2 + 1] + 3 * h[0 + 1] * h[1 + 1]^2,
        3 * h[0 + 1]^2 * h[3 + 1] + 6 * h[0 + 1] * h[1 + 1] * h[2 + 1] + h[1 + 1]^3,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{4}, ::Val{4}) = (
        h[0 + 1]^4,
        4 * h[0 + 1]^3 * h[1 + 1],
        4 * h[0 + 1]^3 * h[2 + 1] + 6 * h[0 + 1]^2 * h[1 + 1]^2,
        4 * h[0 + 1]^3 * h[3 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[2 + 1] +
        4 * h[0 + 1] * h[1 + 1]^3,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{4}, ::Val{5}) = (
        h[0 + 1]^5,
        5 * h[0 + 1]^4 * h[1 + 1],
        5 * h[0 + 1]^4 * h[2 + 1] + 10 * h[0 + 1]^3 * h[1 + 1]^2,
        5 * h[0 + 1]^4 * h[3 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[2 + 1] +
        10 * h[0 + 1]^2 * h[1 + 1]^3,
    )

    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{1}) = (
        h[0 + 1], h[1 + 1], h[2 + 1], h[3 + 1], h[4 + 1]
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{2}) = (
        h[0 + 1]^2,
        2 * h[0 + 1] * h[1 + 1],
        2 * h[0 + 1] * h[2 + 1] + h[1 + 1]^2,
        2 * h[0 + 1] * h[3 + 1] + 2 * h[1 + 1] * h[2 + 1],
        2 * h[0 + 1] * h[4 + 1] + 2 * h[1 + 1] * h[3 + 1] + h[2 + 1]^2,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{3}) = (
        h[0 + 1]^3,
        3 * h[0 + 1]^2 * h[1 + 1],
        3 * h[0 + 1]^2 * h[2 + 1] + 3 * h[0 + 1] * h[1 + 1]^2,
        3 * h[0 + 1]^2 * h[3 + 1] + 6 * h[0 + 1] * h[1 + 1] * h[2 + 1] + h[1 + 1]^3,
        3 * h[0 + 1]^2 * h[4 + 1] +
        6 * h[0 + 1] * h[1 + 1] * h[3 + 1] +
        3 * h[0 + 1] * h[2 + 1]^2 +
        3 * h[1 + 1]^2 * h[2 + 1],
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{4}) = (
        h[0 + 1]^4,
        4 * h[0 + 1]^3 * h[1 + 1],
        4 * h[0 + 1]^3 * h[2 + 1] + 6 * h[0 + 1]^2 * h[1 + 1]^2,
        4 * h[0 + 1]^3 * h[3 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[2 + 1] +
        4 * h[0 + 1] * h[1 + 1]^3,
        4 * h[0 + 1]^3 * h[4 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[3 + 1] +
        6 * h[0 + 1]^2 * h[2 + 1]^2 +
        12 * h[0 + 1] * h[1 + 1]^2 * h[2 + 1] +
        h[1 + 1]^4,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{5}) = (
        h[0 + 1]^5,
        5 * h[0 + 1]^4 * h[1 + 1],
        5 * h[0 + 1]^4 * h[2 + 1] + 10 * h[0 + 1]^3 * h[1 + 1]^2,
        5 * h[0 + 1]^4 * h[3 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[2 + 1] +
        10 * h[0 + 1]^2 * h[1 + 1]^3,
        5 * h[0 + 1]^4 * h[4 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[3 + 1] +
        10 * h[0 + 1]^3 * h[2 + 1]^2 +
        30 * h[0 + 1]^2 * h[1 + 1]^2 * h[2 + 1] +
        5 * h[0 + 1] * h[1 + 1]^4,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{5}, ::Val{6}) = (
        h[0 + 1]^6,
        6 * h[0 + 1]^5 * h[1 + 1],
        6 * h[0 + 1]^5 * h[2 + 1] + 15 * h[0 + 1]^4 * h[1 + 1]^2,
        6 * h[0 + 1]^5 * h[3 + 1] +
        30 * h[0 + 1]^4 * h[1 + 1] * h[2 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1]^3,
        6 * h[0 + 1]^5 * h[4 + 1] +
        30 * h[0 + 1]^4 * h[1 + 1] * h[3 + 1] +
        15 * h[0 + 1]^4 * h[2 + 1]^2 +
        60 * h[0 + 1]^3 * h[1 + 1]^2 * h[2 + 1] +
        15 * h[0 + 1]^2 * h[1 + 1]^4,
    )

    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{1}) = (
        h[0 + 1], h[1 + 1], h[2 + 1], h[3 + 1], h[4 + 1], h[5 + 1]
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{2}) = (
        h[0 + 1]^2,
        2 * h[0 + 1] * h[1 + 1],
        2 * h[0 + 1] * h[2 + 1] + h[1 + 1]^2,
        2 * h[0 + 1] * h[3 + 1] + 2 * h[1 + 1] * h[2 + 1],
        2 * h[0 + 1] * h[4 + 1] + 2 * h[1 + 1] * h[3 + 1] + h[2 + 1]^2,
        2 * h[0 + 1] * h[5 + 1] + 2 * h[1 + 1] * h[4 + 1] + 2 * h[2 + 1] * h[3 + 1],
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{3}) = (
        h[0 + 1]^3,
        3 * h[0 + 1]^2 * h[1 + 1],
        3 * h[0 + 1]^2 * h[2 + 1] + 3 * h[0 + 1] * h[1 + 1]^2,
        3 * h[0 + 1]^2 * h[3 + 1] + 6 * h[0 + 1] * h[1 + 1] * h[2 + 1] + h[1 + 1]^3,
        3 * h[0 + 1]^2 * h[4 + 1] +
        6 * h[0 + 1] * h[1 + 1] * h[3 + 1] +
        3 * h[0 + 1] * h[2 + 1]^2 +
        3 * h[1 + 1]^2 * h[2 + 1],
        3 * h[0 + 1]^2 * h[5 + 1] +
        6 * h[0 + 1] * h[1 + 1] * h[4 + 1] +
        6 * h[0 + 1] * h[2 + 1] * h[3 + 1] +
        3 * h[1 + 1]^2 * h[3 + 1] +
        3 * h[1 + 1] * h[2 + 1]^2,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{4}) = (
        h[0 + 1]^4,
        4 * h[0 + 1]^3 * h[1 + 1],
        4 * h[0 + 1]^3 * h[2 + 1] + 6 * h[0 + 1]^2 * h[1 + 1]^2,
        4 * h[0 + 1]^3 * h[3 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[2 + 1] +
        4 * h[0 + 1] * h[1 + 1]^3,
        4 * h[0 + 1]^3 * h[4 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[3 + 1] +
        6 * h[0 + 1]^2 * h[2 + 1]^2 +
        12 * h[0 + 1] * h[1 + 1]^2 * h[2 + 1] +
        h[1 + 1]^4,
        4 * h[0 + 1]^3 * h[5 + 1] +
        12 * h[0 + 1]^2 * h[1 + 1] * h[4 + 1] +
        12 * h[0 + 1]^2 * h[2 + 1] * h[3 + 1] +
        12 * h[0 + 1] * h[1 + 1]^2 * h[3 + 1] +
        12 * h[0 + 1] * h[1 + 1] * h[2 + 1]^2 +
        4 * h[1 + 1]^3 * h[2 + 1],
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{5}) = (
        h[0 + 1]^5,
        5 * h[0 + 1]^4 * h[1 + 1],
        5 * h[0 + 1]^4 * h[2 + 1] + 10 * h[0 + 1]^3 * h[1 + 1]^2,
        5 * h[0 + 1]^4 * h[3 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[2 + 1] +
        10 * h[0 + 1]^2 * h[1 + 1]^3,
        5 * h[0 + 1]^4 * h[4 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[3 + 1] +
        10 * h[0 + 1]^3 * h[2 + 1]^2 +
        30 * h[0 + 1]^2 * h[1 + 1]^2 * h[2 + 1] +
        5 * h[0 + 1] * h[1 + 1]^4,
        5 * h[0 + 1]^4 * h[5 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1] * h[4 + 1] +
        20 * h[0 + 1]^3 * h[2 + 1] * h[3 + 1] +
        30 * h[0 + 1]^2 * h[1 + 1]^2 * h[3 + 1] +
        30 * h[0 + 1]^2 * h[1 + 1] * h[2 + 1]^2 +
        20 * h[0 + 1] * h[1 + 1]^3 * h[2 + 1] +
        h[1 + 1]^5,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{6}) = (
        h[0 + 1]^6,
        6 * h[0 + 1]^5 * h[1 + 1],
        6 * h[0 + 1]^5 * h[2 + 1] + 15 * h[0 + 1]^4 * h[1 + 1]^2,
        6 * h[0 + 1]^5 * h[3 + 1] +
        30 * h[0 + 1]^4 * h[1 + 1] * h[2 + 1] +
        20 * h[0 + 1]^3 * h[1 + 1]^3,
        6 * h[0 + 1]^5 * h[4 + 1] +
        30 * h[0 + 1]^4 * h[1 + 1] * h[3 + 1] +
        15 * h[0 + 1]^4 * h[2 + 1]^2 +
        60 * h[0 + 1]^3 * h[1 + 1]^2 * h[2 + 1] +
        15 * h[0 + 1]^2 * h[1 + 1]^4,
        6 * h[0 + 1]^5 * h[5 + 1] +
        30 * h[0 + 1]^4 * h[1 + 1] * h[4 + 1] +
        30 * h[0 + 1]^4 * h[2 + 1] * h[3 + 1] +
        60 * h[0 + 1]^3 * h[1 + 1]^2 * h[3 + 1] +
        60 * h[0 + 1]^3 * h[1 + 1] * h[2 + 1]^2 +
        60 * h[0 + 1]^2 * h[1 + 1]^3 * h[2 + 1] +
        6 * h[0 + 1] * h[1 + 1]^5,
    )
    test_hⁱ_mod_xⁿ⁻¹(h::NTuple{6}, ::Val{7}) = (
        h[0 + 1]^7,
        7 * h[0 + 1]^6 * h[1 + 1],
        7 * h[0 + 1]^6 * h[2 + 1] + 21 * h[0 + 1]^5 * h[1 + 1]^2,
        7 * h[0 + 1]^6 * h[3 + 1] +
        42 * h[0 + 1]^5 * h[1 + 1] * h[2 + 1] +
        35 * h[0 + 1]^4 * h[1 + 1]^3,
        7 * h[0 + 1]^6 * h[4 + 1] +
        42 * h[0 + 1]^5 * h[1 + 1] * h[3 + 1] +
        21 * h[0 + 1]^5 * h[2 + 1]^2 +
        105 * h[0 + 1]^4 * h[1 + 1]^2 * h[2 + 1] +
        35 * h[0 + 1]^3 * h[1 + 1]^4,
        7 * h[0 + 1]^6 * h[5 + 1] +
        42 * h[0 + 1]^5 * h[1 + 1] * h[4 + 1] +
        42 * h[0 + 1]^5 * h[2 + 1] * h[3 + 1] +
        105 * h[0 + 1]^4 * h[1 + 1]^2 * h[3 + 1] +
        105 * h[0 + 1]^4 * h[1 + 1] * h[2 + 1]^2 +
        140 * h[0 + 1]^3 * h[1 + 1]^3 * h[2 + 1] +
        21 * h[0 + 1]^2 * h[1 + 1]^5,
    )

    for N ∈ 1:6
        h = (3, 7, 11, 17, 29, 37)[1:N]
        hⁱ⁻¹ = (1, 0, 0, 0, 0, 0)[1:N]
        for i ∈ 1:(N + 1)
            hⁱ = hⁱ✖h_mod_xⁿ⁻¹(hⁱ⁻¹, h)
            @test all(hⁱ .≈ test_hⁱ_mod_xⁿ⁻¹(h, Val(i)))
            hⁱ⁻¹ = hⁱ
        end
    end
end

@doc raw"""
    lagrange_inversion(a)

Compute the *compositional* inverse (a.k.a. reversion) of the power series expansion
```math
f(x) = a[1]*x + a[2]*x^2 + \ldots + a[n-1]*x^{n-1}
```
about 0, where `a` is an NTuple.  Note, in particular, that there is no constant term.  The
result is a similar NTuple `b` allowing us to write
```math
f^{-1}(y) = b[1]*y + b[2]*y^2 + \ldots + b[n-1]*y^{n-1},
```
such that ``f^{-1}(f(x)) = f(f^{-1}(y)) = x`` mod ``x^n``.

See [`truncated_series_inverse`](@ref) for the *multiplicative* inverse.

When the constant coefficient ``a_0`` is nonzero, the result must be expanded about a
different point, which is done by evaluating the output as ``f^{-1}(y-a_0)``.  Similarly, if
the original expansion is about a point ``x_0 ≠ 0``, the result must be shifted by adding
``x_0`` to the output.

[Johansson (2015)](https://doi.org/10.1090/S0025-5718-2014-02857-3) summarizes this basic
form of the algorithm nicely:

> Our setting is the ring of truncated power series ``R[[x]]/\langle x^n \rangle`` over a
> commutative coefficient ring ``R`` in which the integers ``1,...,n−1`` are cancellable
> (i.e., nonzero and not zero divisors).  [...]  If ``f(x) = x/h(x)`` where ``h(0)`` is a
> unit, then the compositional inverse or reversion ``f^{-1}(x)`` satisfying ``f(f^{-1}(x))
> = f^{-1}(f(x)) = x`` exists and its coefficients are given by
> ```math
> [x^k] f^{-1}(x) = \frac{1}{k} [x^{k-1}] h(x)^k.
> ```

Note that Johansson also presents a pair of asymptotically faster algorithms for computing
the compositional inverse.  Because of the low orders of the power series expansions we
typically work with, it is not clear if the improved scaling of those algorithms would
actually be beneficial in practice, so we stick with the basic algorithm here — though they
would not be too difficult to implement if needed.

"""
function lagrange_inversion(a::NTuple{N,T}) where {N,T}
    h = x╱f_mod_xⁿ⁻¹(a)
    f⁻¹ = zeros(MVector{N,typeof(h[end] / 2)})
    hⁱ = h  # Create storage for the loop
    for i ∈ eachindex(f⁻¹)
        f⁻¹[i] = hⁱ[i] / i  # Note that hⁱ[i] is the coefficient of x^(i-1) in hⁱ(x)
        hⁱ = hⁱ✖h_mod_xⁿ⁻¹(hⁱ, h)  # We set hⁱ to hⁱ⁺¹ for the next iteration
    end
    return Tuple(f⁻¹)
end

@testitem "lagrange_inversion" begin
    using PostNewtonian: lagrange_inversion
    using Random

    rng = Random.Xoshiro(123)
    for N ∈ 2:16
        n = N + 1
        a = Tuple(randn(rng, N))
        b = lagrange_inversion(a)
        for x ∈ (1e-12, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 8e-4, 4e-4, 2e-4, 1e-4)
            f = evalpoly(x, (zero(eltype(a)), a...))
            g = evalpoly(f, (zero(eltype(b)), b...))
            δ = abs(x - g)
            ϵ = max(n * eps(eltype(a)), 4n * x^n)
            # println()
            # @show n N x g x^n δ ϵ
            @test δ < ϵ
        end
        # println()
    end
end
