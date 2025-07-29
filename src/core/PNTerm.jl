# We enclose all of this in a `baremodule` so that we can isolate Base operations, and use
# PNBase operations by default.
baremodule PNTerms

# See explanation in `PNBase.jl` for why we use these operators.
using Base: Base, @doc, @assert, one, zero, >, ÷, + as ⊕, - as ⊖, * as ⊛, / as ⊘, ^ as ↑

import PostNewtonian: constant_convert
using PostNewtonian: PostNewtonian, PNSystem
using PostNewtonian.InlineExports: @export
using PostNewtonian.PNBase


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
@export struct PNTerm{T,PNOrder,c⁻¹Exponent}
    coeff::T

    function PNTerm{T,PNOrder,c⁻¹Exponent}(coeff) where {T,PNOrder,c⁻¹Exponent}
        if c⁻¹Exponent > 2 ⊛ PNOrder
            coeff = zero(coeff)
        end
        return new{T,PNOrder,c⁻¹Exponent}(coeff)
    end
    function PNTerm{T,PNOrder}(c⁻¹exp::Int, coeff) where {T,PNOrder}
        if c⁻¹exp > 2 ⊛ PNOrder
            coeff = zero(coeff)
        end
        return new{T,PNOrder,c⁻¹exp}(coeff)
    end
end

"""
    PNExpansionParameter(pnsystem)

Create a [`PNTerm`](@ref) object representing the post-Newtonian expansion parameter ``c``.
This can be used to automatically create more complicated `PNTerm`s, which combine to form a
[`PNExpansion`](@ref).  This is a simple but effective way to write PN formulas while
automatically tracking the PN order of each term.
"""
@export function PNExpansionParameter(::T) where {NT,PNOrder,T<:PNSystem{NT,PNOrder}}
    return PNTerm{NT,PNOrder}(-1, one(T))
end

Base.length(pn::PNTerm) = 1  # COV_EXCL_LINE  (Covered, but missed due to const-folding)
Base.eltype(pn::PNTerm{T}) where {T} = T

@export c⁻¹exp(pn::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent} =
    c⁻¹Exponent

function PostNewtonian.constant_convert(
    pn::P, term::T
) where {NT,PNOrder,P<:PNSystem{NT,PNOrder},T<:PNTerm{NT,PNOrder}}
    term
end

### NOTE: Adding or subtracting multiple `PNTerm`s will produce a `PNExpansion`, which
### is defined in `PNExpansion.jl`.  The following additive functions are just unary.

function Base.sum(term::PNTerm)
    return term.coeff
end

function PNBase.:+(term::PNTerm)
    return term
end

function PNBase.:-(term::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent}
    return PNTerm{T,PNOrder,c⁻¹Exponent}(⊖(term.coeff))
end

function Base.inv(term::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent}
    return PNTerm{T,PNOrder}(⊖(c⁻¹Exponent), Base.inv(term.coeff))
end

function PNBase.:√(term::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent}
    @assert Base.iseven(c⁻¹Exponent) "Only half-integer PN orders are supported."
    return PNTerm{T,PNOrder}(c⁻¹Exponent ÷ 2, Base.sqrt(term.coeff))
end

function PNBase.:^(
    term::PNTerm{T,PNOrder,c⁻¹Exponent}, n::Int
) where {T,PNOrder,c⁻¹Exponent}
    coeff = term.coeff↑n
    return PNTerm{typeof(coeff),PNOrder}(c⁻¹exp(term) ⊛ n, coeff)
end

function PNBase.:*(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    c⁻¹Exponent = c⁻¹exp(term1) ⊕ c⁻¹exp(term2)
    coeff = term1.coeff ⊛ term2.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

function PNBase.:*(
    x::Number, term::PNTerm{T,PNOrder,c⁻¹Exponent}
) where {T,PNOrder,c⁻¹Exponent}
    coeff = x ⊛ term.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

PNBase.:*(term::PNTerm, x::Number) = PNBase.:*(x, term)

function PNBase.:/(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    c⁻¹Exponent = c⁻¹E1 ⊖ c⁻¹E2
    coeff = term1.coeff ⊘ term2.coeff
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

function PNBase.:/(
    term::PNTerm{T,PNOrder,c⁻¹Exponent}, x::Number
) where {T,PNOrder,c⁻¹Exponent}
    coeff = term.coeff ⊘ x
    return PNTerm{typeof(coeff),PNOrder,c⁻¹Exponent}(coeff)
end

function PNBase.:/(
    x::Number, term::PNTerm{T,PNOrder,c⁻¹Exponent}
) where {T,PNOrder,c⁻¹Exponent}
    coeff = x ⊘ term.coeff
    return PNTerm{typeof(coeff),PNOrder}(⊖(c⁻¹exp(term)), coeff)
end

end  # baremodule PNTerms

@testitem "PNTerm algebra" begin
    using Base: Base, one, zero, <, ÷, + as ⊕, - as ⊖, * as ⊛, / as ⊘, ^ as ↑
    using DoubleFloats: Double64
    using PostNewtonian: PostNewtonian, pn_order
    using PostNewtonian.PNBase: ln, (√), (+), (-), (*), (/), (//), (^)
    using PostNewtonian.PNTerms: PNTerm, PNExpansionParameter, c⁻¹exp, PNBase, constant_convert
    using PostNewtonian.PNExpansions: PNExpansion

    for T ∈ [Float64, Float16, Double64]
        pn = BBH(randn(T, 14), 9//2)
        pn[:v] = /(pn, 1, 5)
        c = PNExpansionParameter(pn)
        z = zero(T)
        x = /(pn, 6, 5)
        y = /(pn, 17, 5)
        w = /(pn, 28, 5)
        v = /(pn, 27, 13)

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
            @test c⁻¹exp(term) == c⁻¹exponent
            @test term.coeff == coeff
            @test term.coeff isa eltype(pn)
            @test length(term) == 1
            @test eltype(term) ≡ T
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
            @test eltype(expr) ≡ T
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

    pn = BBH(randn(14), 2//1)

    # constant_convert should just return the term
    t0 = PNTerm{Float64,2//1}(0, 5.0)
    @test constant_convert(pn, t0) ≡ t0

    # Base.sum
    @test sum(t0) == 5.0

    # unary + (PNBase.+)
    t1 = PNTerm{Float64,2}(1, 2.5)
    @test +t1 ≡ t1

    # unary - (PNBase.-)
    t1n = -t1
    @test c⁻¹exp(t1n) == c⁻¹exp(t1)
    @test t1n.coeff == -2.5

    # inv
    ti = inv(t1)
    @test c⁻¹exp(ti) == -c⁻¹exp(t1)
    @test ti.coeff == 1/2.5

    # sqrt on even exponent
    t2 = PNTerm{Float64,2}(2, 4.0)
    ts = √(t2)
    @test c⁻¹exp(ts) == 1
    @test ts.coeff == 2.0

    # sqrt on odd exponent should error
    t3 = PNTerm{Float64,2}(1, 9.0)
    @test_throws AssertionError √(t3)

    # power
    tp = t1 ^ 3
    @test c⁻¹exp(tp) == 3*c⁻¹exp(t1)
    @test tp.coeff == 2.5^3

    # multiply term * term
    ta = PNTerm{Float64,2}(1, 2.0)
    tb = PNTerm{Float64,2}(2, 3.0)
    tab = ta * tb
    @test c⁻¹exp(tab) == 3
    @test tab.coeff == 6.0

    # scalar * term and term * scalar
    ts1 = 4.0 * ta
    ts2 = ta * 4.0
    @test c⁻¹exp(ts1) == c⁻¹exp(ta)
    @test ts1.coeff == 8.0
    @test ts2.coeff == 8.0

    # term / term
    td = tb / ta
    @test c⁻¹exp(td) == 1      # 2 - 1
    @test td.coeff == 1.5

    # term / scalar
    td1 = tb / 2.0
    @test c⁻¹exp(td1) == c⁻¹exp(tb)
    @test td1.coeff == 1.5

    # scalar / term
    td2 = 12.0 / ta
    @test c⁻¹exp(td2) == -1    # -c⁻¹exp(ta)
    @test td2.coeff == 6.0
end
