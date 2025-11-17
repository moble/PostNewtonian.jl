# We enclose all of this in a `baremodule` so that we can isolate Base operations, and use
# PNBase operations by default.
baremodule PNExpansions

# See explanation in `PNBase.jl` for why we use these operators.
using Base: Base, NTuple, ntuple, Val, Vector, @doc, @raw_str, @assert, @inbounds, @__dot__,
    (:), eltype, promote_type, isbitstype, min, max, one, zero,
    ==, ≤, <, >, ÷, + as ⊕, - as ⊖, * as ⊛, / as ⊘, ^ as ↑

using PostNewtonian: PostNewtonian
using PostNewtonian.PNBase
using PostNewtonian.PNTerms
using PostNewtonian.InlineExports: @export, @public
using StaticArrays: MVector, SVector
import FastDifferentiation

# This a utility that allow us to interoperate with FastDifferentiation.Node and other
# Number types.
function _efficient_vector(::Val{N}, ::Val{T}) where {N,T}
    if isbitstype(T)
        MVector{N,T}(undef)
    else
        Vector{T}(undef, N)
    end
end

@doc raw"""
    PNExpansion{N,T,NMax}

This object can be multiplied by a scalar or another `PNExpansion` object, and contains a
tuple of coefficients.  The coefficients are stored in the order of the expansion, with the
zeroth-order coefficient first.  Multiplication by a scalar just multiplies each of the
elements.  Multiplication by another `PNExpansion` object is done by convolving the two
tuples of coefficients.

Blanchet (2014) defines the post-Newtonian expansion parameter as follows:

> This parameter represents essentially a slow motion estimate ``ϵ ∼ 𝑣/𝑐``, where ``𝑣``
> denotes a typical internal velocity.  By a slight abuse of notation, following
> Chandrasekhar et al. [...], we shall henceforth write formally ``ϵ ≡ 1/𝑐``, even though
> ``ϵ`` is dimensionless whereas ``𝑐`` has the dimension of a velocity. Thus, ``1/𝑐 ≪ 1``
> in the case of post-Newtonian sources. The small post-Newtonian remainders will be denoted
> ``𝒪(1/𝑐^𝑛)``. Furthermore, [...] we shall refer to a small post-Newtonian term with
> formal order ``𝒪(1/𝑐^𝑛)`` relative to the Newtonian acceleration in the equations of
> motion, as ``\frac{𝑛}{2}\text{PN}``.

Therefore, we consider the coefficients of the `PNExpansion` to be a polynomial in ``1/𝑐``.
Here, the type parameter `N` corresponds to the number of elements actually present in the
tuple of coefficients, and `T` is the type of the coefficients.  The `NMax` parameter is the
maximum number of elements, related to the usual PN order by
```math
\text{pn_order} = \frac{\texttt{NMax}-1} {2}.
```
The `N` parameter is not related to the PN order; it is just used by Julia to know how many
elements are currently in the coefficients, but is required to be 1 ≤ N ≤ NMax.

"""
@public struct PNExpansion{N,T,NMax}
    coeffs::NTuple{N,T}

    function PNExpansion{N,T,NMax}(coeffs) where {N,T,NMax}
        if N < 1
            throw(ArgumentError("`N=$N` must be greater than 0."))  # COV_EXCL_LINE
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be less than `NMax=$NMax`."))  # COV_EXCL_LINE
        end
        return new{N,T,NMax}(coeffs)
    end
    function PNExpansion(coeffs::NTuple{N,T}, NMax) where {N,T}
        if N < 1
            throw(ArgumentError("`N=$N` must be greater than 0."))
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be less than `NMax=$NMax`."))  # COV_EXCL_LINE
        end
        return new{N,T,NMax}(coeffs)
    end
end

Base.Tuple(pn::PNExpansion) = pn.coeffs
SVector(pn::PNExpansion) = SVector(pn.coeffs)

PostNewtonian.pn_order(::PNExpansion{N,T,NMax}) where {N,T,NMax} = (NMax ⊖ 1)//2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = Base.length(pn.coeffs)
Base.eltype(pn::PNExpansion) = Base.eltype(pn.coeffs)

function Base.sum(pn_expansion::PNExpansion{N,T,NMax}) where {N,T,NMax}
    return Base.sum(pn_expansion.coeffs; init=zero(T))
end

function PNBase.:+(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    return PNExpansion(ntuple(i -> i == 1 ? pn[1] + x : T3(pn[i]), Val(N)), NMax)
end

PNBase.:+(x::T, pn::PNExpansion) where {T<:Number} = PNBase.:+(pn, x)

function PNBase.:+(
    pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}
) where {N1,N2,T1,T2,NMax1,NMax2}
    throw(
        ArgumentError(
            "`PNExpansion` addition is only defined for objects of the same PN order." ⊛
            "\nGot NMax1=$(NMax1) and NMax2=$(NMax2).",
        ),
    )
end

function PNBase.:+(
    pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    if N1 > N2
        return pn2 + pn1
    else
        PNExpansion(ntuple(i -> sum_term(i, pn1, pn2), Val(N2)), NMax)
    end
end

function sum_term(
    i, pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    T3 = promote_type(T1, T2)
    if i ≤ N1
        return T3(pn1.coeffs[i] ⊕ pn2.coeffs[i])
    else
        return T3(pn2.coeffs[i])
    end
end

function PNBase.:+(
    x::T1, term::PNTerm{T2,PNOrder,c⁻¹Exponent}
) where {T1<:Number,T2,PNOrder,c⁻¹Exponent}
    if c⁻¹exp(term) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " ⊛
                "c⁻¹exp(term)=$(c⁻¹exp(term))." ⊛
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    T = promote_type(T1, T2)
    N₀ = c⁻¹exp(term) ⊕ 1
    NMax = Int(2 ⊛ PNOrder ⊕ 1)
    N = min(N₀, NMax)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds coeffs[1] = x
    @inbounds if N₀ ≤ NMax
        coeffs[N₀] += term.coeff
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end

PNBase.:+(term::PNTerm, x::Number) = PNBase.:+(x, term)

function PNBase.:+(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    if c⁻¹exp(term1) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " ⊛
                "c⁻¹exp(term1)=$(c⁻¹exp(term1))." ⊛
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    if c⁻¹exp(term2) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " ⊛
                "c⁻¹exp(term2)=$(c⁻¹exp(term2))." ⊛
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    T = promote_type(T1, T2)
    N1₀ = c⁻¹exp(term1) ⊕ 1
    N2₀ = c⁻¹exp(term2) ⊕ 1
    NMax = Int(2 ⊛ PNOrder ⊕ 1)
    N = min(max(N1₀, N2₀), NMax)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds if N1₀ ≤ N
        coeffs[N1₀] += term1.coeff
    end
    @inbounds if N2₀ ≤ N
        coeffs[N2₀] += term2.coeff
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end

function PNBase.:+(
    term::PNTerm{T1,PNOrder,c⁻¹E1}, expansion::PNExpansion{N2,T2,NMax2}
) where {T1,PNOrder,c⁻¹E1,N2,T2,NMax2}
    if c⁻¹exp(term) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " ⊛
                "c⁻¹exp(term)=$(c⁻¹exp(term))." ⊛
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    N1 = c⁻¹exp(term) ⊕ 1
    NMax1 = Int(2 ⊛ PNOrder ⊕ 1)
    NMax = min(NMax1, NMax2)
    N = min(max(N1, N2), NMax)
    T = promote_type(T1, T2)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds if N1 ≤ N
        coeffs[N1] += term.coeff
    end
    @inbounds for i ∈ 1:N
        if i ≤ N2
            coeffs[i] += expansion[i]
        end
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end

PNBase.:+(expansion::PNExpansion, term::PNTerm) = PNBase.:+(term, expansion)

PNBase.:-(term1::PNTerm, term2::PNTerm) = PNBase.:+(term1, -term2)
PNBase.:-(term::PNTerm, x::Number) = PNBase.:+(term, -x)
PNBase.:-(x::Number, term::PNTerm) = PNBase.:+(x, -term)
PNBase.:-(term::PNTerm, expansion::PNExpansion) = PNBase.:+(term, -expansion)
PNBase.:-(expansion::PNExpansion, term::PNTerm) = PNBase.:+(expansion, -term)
PNBase.:-(x::Number, expansion::PNExpansion) = PNBase.:+(x, -expansion)
PNBase.:-(expansion::PNExpansion, x::Number) = PNBase.:+(expansion, -x)

function PNBase.:-(pn::PNExpansion{N,T,NMax}) where {N,T,NMax}
    return PNExpansion{N,T,NMax}((-).(pn.coeffs))
end

function PNBase.:*(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    return PNExpansion{N,T3,NMax}(@. T3(pn.coeffs * x))
end

PNBase.:*(x::T, pn::PNExpansion) where {T<:Number} = PNBase.:*(pn, x)

function PNBase.:*(
    pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}
) where {N1,N2,T1,T2,NMax1,NMax2}
    throw(
        ArgumentError(
            "`PNExpansion` multiplication is only defined for objects of the same PN order." ⊛
            "\nGot NMax1=$(NMax1) and NMax2=$(NMax2).",
        ),
    )
end

function PNBase.:*(
    pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    if N1 > N2
        return PNBase.:*(pn2, pn1)
    else
        N3 = min(N1 + N2 - 1, NMax)
        PNExpansion(ntuple(i -> product_term(i, pn1, pn2), Val(N3)), NMax)
    end
end

function product_term(
    i, pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    T3 = promote_type(T1, T2)
    return Base.sum(
        pn1.coeffs[j] ⊛ pn2.coeffs[i - j + 1] for j ∈ max(1, i - N2 + 1):min(i, N1);
        init=zero(T3),
    )
end

function PNBase.:*(
    expansion::PNExpansion{N1,T1,NMax1}, term::PNTerm{T2,PNOrder,c⁻¹E2}
) where {N1,T1,NMax1,T2,PNOrder,c⁻¹E2}
    ΔN = c⁻¹exp(term)  # Note that ΔN may be negative!
    NMax2 = Int(2 ⊛ PNOrder ⊕ 1)
    NMax = min(NMax1, NMax2)
    N = min(max(N1, N1 ⊕ ΔN), NMax)

    # Check that no terms from expansion will be lost to negative PN orders
    @inbounds for i ∈ 1:min(max(0, ⊖(ΔN)), N1)
        if !iszero(expansion[i])
            throw(
                ArgumentError(
                    "Cannot multiply `PNExpansion` by `PNTerm` with negative exponent: " ⊛
                    "c⁻¹exp(term)=$(c⁻¹exp(term)).\nThe result will be a `PNExpansion`, " ⊛
                    "which cannot store positive exponents.",
                ),
            )
        end
    end

    T = promote_type(T1, T2)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds for i ∈ max(1, 1 ⊖ ΔN):min(N1, N ⊖ ΔN)
        coeffs[i ⊕ ΔN] = expansion[i] ⊛ term.coeff
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end

PNBase.:*(term::PNTerm, expansion::PNExpansion) = PNBase.:*(expansion, term)
# (a, b, c, d, e, f, g) * (c⁻¹^2) = (0, 0, a, b, c, d, e)

PNBase.:/(expansion::PNExpansion, x) = PNBase.:*(expansion, Base.inv(x))

# Note that an PNExpansion is really a *relative* expansion in `1/c` — *not* `v` or `x`.
# Therefore, the correct derivative with respect to any variable (other than `c`, which we
# never differentiate with respect to) extends to just derivatives of the coefficients,
# without any change to the exponent of `1/c`.
function FastDifferentiation.derivative(
    pn_expansion::PNExpansion{N,T,NMax}, fd_node::FastDifferentiation.Node
) where {N,T,NMax}
    return PNExpansion(
        ntuple(i -> FastDifferentiation.derivative(pn_expansion[i], fd_node), Val(N)), NMax
    )
end

end  # baremodule PNExpansions

@testitem "PNExpansion algebra" begin
    using Base: Base, one, zero, <, ÷, + as ⊕, - as ⊖, * as ⊛, / as ⊘, ^ as ↑
    using StaticArrays: MVector, SVector
    using Symbolics: @variables, simplify, substitute
    using PostNewtonian.PNBase: ln, (√), (+), (-), (*), (/), (//), (^)
    using PostNewtonian.PNExpansions: PNExpansion

    # Test edge cases
    @test_throws ArgumentError PNExpansion{0, Float64, 1}(())
    @test_throws ArgumentError PNExpansion{2, Float64, 1}((1.2, 3.4))
    @test_throws ArgumentError PNExpansion((), 0)
    @test_throws ArgumentError PNExpansion((1.2, 3.4), 1)

    for N1 ∈ 1:9
        # Test conversions
        coeffs = ntuple(i -> i + 1.0, N1)
        for NMax ∈ N1:(N1 ⊕ 3)
            expansion = PNExpansion(coeffs, NMax)
            @test Base.Tuple(expansion) == coeffs
            @test SVector(expansion) == SVector(coeffs)
            @test Base.eltype(expansion) == Float64
            @test Base.length(expansion) == N1
            for i ∈ 1:N1
                @test expansion[i] == coeffs[i]
            end
        end

        for N2 ∈ 1:9
            for NMax ∈ max(N1, N2):(N1 ⊕ N2 ⊕ 3)

                @variables c⁻¹ x[1:N1] y[1:N2] z
                poly(e::PNExpansion) = sum(e[i] ⊛ c⁻¹↑(i ⊖ 1) for i ∈ 1:length(e))
                eˣ = PNExpansion(tuple(x...), NMax)
                eʸ = PNExpansion(tuple(y...), NMax)

                # Test sums
                polysum = simplify(poly(eˣ + eʸ); expand=true)
                sumpoly = simplify(poly(eˣ) ⊕ poly(eʸ); expand=true)
                Δ = simplify(polysum ⊖ sumpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ + PNExpansion(tuple(z, x...), NMax ⊕ 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax ⊕ 1) + eˣ

                # Test products
                polyprod = simplify(poly(eˣ * eʸ); expand=true)
                prodpoly = simplify(
                    substitute(
                        simplify(poly(eˣ) ⊛ poly(eʸ); expand=true),
                        Dict([c⁻¹↑n => 0 for n ∈ NMax:(2NMax ⊕ 3)]),
                    );
                    expand=true,
                )
                Δ = simplify(polyprod ⊖ prodpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ * PNExpansion(tuple(z, x...), NMax ⊕ 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax ⊕ 1) * eˣ
            end
        end
    end
end
