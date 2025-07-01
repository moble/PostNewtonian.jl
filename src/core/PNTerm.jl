"""
    PNExpansionParameter(pnsystem)

Create a [`PNTerm`](@ref) object representing the post-Newtonian expansion parameter ``c``.
This can be used to automatically create more complicated `PNTerm`s, which combine to form a
[`PNExpansion`](@ref).  This is a simple but effective way to write PN formulas while
automatically tracking the PN order of each term.
"""
@public function PNExpansionParameter(::T) where {NT,PNOrder,T<:PNSystem{NT,PNOrder}}
    return PNTerm{NT,PNOrder}(-1, one(T))
end

"""
    PNTerm{T,PNOrder,c⁻¹Exponent}

This object represents a single term in a PNExpansion.  It has a single field: `coeff`,
which is the coefficient of the term.  The type parameter `T` is the type of the
coefficient.  The type parameter `PNOrder` is a half-integer (just as in
[`PNSystem`](@ref)s) representing the PN order of the expansion.  And the type parameter
`c⁻¹Exponent` is an integer representing the exponent of the PN expansion parameter ``1/c``.

`PNTerm`s can be multiplied and divided by scalars and exponentiated by integers, to produce
another `PNTerm`.  They can also be added to other `PNTerm`s to produce a `PNExpansion`.

A simple way to define a `PNTerm` or a `PNExpansion` is to define the PN expansion parameter
```julia
c = PNExpansionParameter(pnsystem)
```
and use that naturally in formulas, as in
```julia
e = 1 + (v/c)^2 * (-ν/12 - 3//4) + (v/c)^4 * (-ν^2/24 + 19ν/8 - 27//8)
```
Any exponent higher than the desired `PNOrder` will be automatically set to zero.

Useful facts:
  - `v` has order `1/c`
  - `x` has order `1/c^2`
  - `γ` has order `1/c^2`
  - `1/r` has order `1/c^2`

"""
@public struct PNTerm{T,PNOrder,c⁻¹Exponent}
    coeff::T

    function PNTerm{T,PNOrder,c⁻¹Exponent}(coeff) where {T,PNOrder,c⁻¹Exponent}
        if c⁻¹Exponent > 2PNOrder
            coeff = zero(coeff)
        end
        return new{T,PNOrder,c⁻¹Exponent}(coeff)
    end
    function PNTerm{T,PNOrder}(c⁻¹exp::Int, coeff) where {T,PNOrder}
        if c⁻¹exp > 2PNOrder
            coeff = zero(coeff)
        end
        return new{T,PNOrder,c⁻¹exp}(coeff)
    end
end

Base.length(pn::PNTerm) = 1
Base.eltype(pn::PNTerm{T}) where {T} = T
@public c⁻¹exp(pn::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent} =
    c⁻¹Exponent

function Base.sum(pn::PNTerm)
    return pn.coeff
end

function Base.:+(pn::PNTerm)
    return pn
end

function Base.inv(term::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent}
    return PNTerm{T,PNOrder}(-c⁻¹exp(term), inv(term.coeff))
end

function Base.:^(term::PNTerm{T,PNOrder,c⁻¹Exponent}, n::Int) where {T,PNOrder,c⁻¹Exponent}
    coeff = term.coeff^n
    return PNTerm{typeof(coeff),PNOrder}(c⁻¹exp(term) * n, coeff)
end

function Base.:*(
    x::Number, term::PNTerm{T,PNOrder,c⁻¹Exponent}
) where {T,PNOrder,c⁻¹Exponent}
    coeff = x * term.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end
Base.:*(term::PNTerm, x::Number) = x * term

function Base.:/(
    term::PNTerm{T,PNOrder,c⁻¹Exponent}, x::Number
) where {T,PNOrder,c⁻¹Exponent}
    coeff = term.coeff / x
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

function Base.:/(
    x::Number, term::PNTerm{T,PNOrder,c⁻¹Exponent}
) where {T,PNOrder,c⁻¹Exponent}
    coeff = x / term.coeff
    return PNTerm{typeof(coeff),PNOrder}(-c⁻¹exp(term), coeff)
end

function Base.:*(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    c⁻¹Exponent = c⁻¹exp(term1) + c⁻¹exp(term2)
    coeff = term1.coeff * term2.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

function Base.:/(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    c⁻¹Exponent = c⁻¹E1 - c⁻¹E2
    coeff = term1.coeff / term2.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

@testitem "PNTerm algebra" begin
    using DoubleFloats
    using PostNewtonian: PostNewtonian, PNExpansionParameter, PNExpansion, pn_order

    for T ∈ [Float64, Float16, Double64]
        pn = BBH(randn(T, 14), 9//2)
        pn[:v] = T(1) / 5
        c = PNExpansionParameter(pn)
        z = zero(T)
        x = T(6) / 5
        y = T(17) / 5
        w = T(28) / 5
        v = T(27) / 13

        # Test behavior of `c` as the basic PNTerm
        for (term, c⁻¹exponent, coeff) ∈ (
            (c, -1, 1),
            (c^2, -2, 1),
            (x * c^2, -2, x),
            (c^2 * x, -2, x),
            (c^2 / x, -2, 1 / x),
            ((x * c)^2, -2, x^2),
            ((x * c^2) / c^4, 2, x),
            ((x * c^2) / c^-4, -6, x),
            ((x / c^2) / c^4, 6, x),
            ((x / c^2) * (y / c^4), 6, x * y),
        )
            @test PostNewtonian.c⁻¹exp(term) == c⁻¹exponent
            @test term.coeff == coeff
            @test term.coeff isa eltype(pn)
            @test length(term) == 1
            @test eltype(term) === T
        end

        # Test PNExpressions
        @test_throws ArgumentError PNExpansion((), 0)
        for (expr, expected) ∈ (
            (w - (x - y / c), (w - x, y)),
            (w - (x / c - y), (w + y, -x)),
            (x - y / c, (x, -y)),
            (x / c - y, (-y, x)),
            (x / c^6 + y / c, (z, y, z, z, z, z, x)),
            (x / c^6 + y / c + w / c^10, (z, y, z, z, z, z, x, z, z, z)),
            (x / c^6 - y / c, (z, -y, z, z, z, z, x)),
            (x / c^6 - y / c - w, (-w, -y, z, z, z, z, x)),
            (x / c^6 - y / c - w / c^10, (z, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - (y / c - w), (w, -y, z, z, z, z, x)),
            (x / c^6 - (y / c - w / c^10), (z, -y, z, z, z, z, x, z, z, z)),
            (-(x / c^6) - y / c, (z, -y, z, z, z, z, -x)),
            (-(x / c^6) - y / c - w, (-w, -y, z, z, z, z, -x)),
            (-(x / c^6) - y / c - w / c^10, (z, -y, z, z, z, z, -x, z, z, z)),
            ((x * c^2) / c^4 + y / c, (z, y, x)),
            ((x * c^2) / c^4 + y / c + w, (w, y, x)),
            ((x * c^2) / c^9 + y / c + w, (w, y, z, z, z, z, z, x)),
            (w + (x * c^2) / c^4 + y / c, (w, y, x)),
            (w + (x * c^2) / c^9 + y / c, (w, y, z, z, z, z, z, x)),
            (w + ((x * c^2) / c^4 + y / c), (w, y, x)),
            (w + ((x * c^2) / c^9 + y / c), (w, y, z, z, z, z, z, x)),
            (((x * c^2) / c^4 + y / c + w) / c^3, (z, z, z, w, y, x)),
            (((x * c^2) / c^4 + y / c + w) / c^5, (z, z, z, z, z, w, y, x)),
            (((x * c^2) / c^7 + y / c + w) / c^5, (z, z, z, z, z, w, y, z, z, z)),
            (v * (((x * c^2) / c^4 + y / c + w) / c^3), v .* (z, z, z, w, y, x)),
            ((((x * c^2) / c^4 + y / c + w) / c^3) * v, v .* (z, z, z, w, y, x)),
            ((((x * c^2) / c^4 + y / c + w) / c^3) / v, (z, z, z, w, y, x) .* (1 / v)),
        )
            @test expr.coeffs == expected
            @test pn_order(expr) == pn_order(pn)
            @test sum(expr) == sum(expected)
            @test eltype(expr) === T
        end

        # Can't make a PNExpression with positive exponents
        @test_throws ArgumentError x * c + y
        @test_throws ArgumentError x * c + y / c
        @test_throws ArgumentError x * c^2 + y
        @test_throws ArgumentError x * c^2 + y / c
        @test_throws ArgumentError y + x * c
        @test_throws ArgumentError y / c + x * c
        @test_throws ArgumentError y + x * c^2
        @test_throws ArgumentError y / c + x * c^2
        @test_throws ArgumentError x * c^2 + (y / c + z / c^2)
    end
end
