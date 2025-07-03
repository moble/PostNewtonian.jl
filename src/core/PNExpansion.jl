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
            throw(ArgumentError("`N=$N` must be >0."))
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be <`NMax=$NMax`."))
        end
        return new{N,T,NMax}(coeffs)
    end
    function PNExpansion(coeffs::NTuple{N,T}, NMax) where {N,T}
        if N < 1
            throw(ArgumentError("`N=$N` must be >0."))
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be <`NMax=$NMax`."))
        end
        return new{N,T,NMax}(coeffs)
    end
end

pn_order(::PNExpansion{N,T,NMax}) where {N,T,NMax} = (NMax - 1)//2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = length(pn.coeffs)
Base.eltype(pn::PNExpansion) = eltype(pn.coeffs)

function Base.sum(pn_expansion::PNExpansion{N,T,NMax}) where {N,T,NMax}
    return sum(pn_expansion[i] for i ∈ 1:N, init ∈ zero(T))
end

function Base.:+(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    return PNExpansion(ntuple(i -> i == 1 ? pn[1] + x : T3(pn[i]), Val(N)), NMax)
end
Base.:+(x::T, pn::PNExpansion) where {T<:Number} = pn + x

function Base.:-(pn::PNExpansion{N,T,NMax}) where {N,T,NMax}
    return PNExpansion{N,T,NMax}((-).(pn.coeffs))
end

function Base.:*(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    return PNExpansion{N,T3,NMax}(@. T3(pn.coeffs * x))
end
Base.:*(x::T, pn::PNExpansion) where {T<:Number} = pn * x

function Base.:+(
    pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}
) where {N1,N2,T1,T2,NMax1,NMax2}
    throw(
        ArgumentError(
            "`PNExpansion` addition is only defined for objects of the same PN order." *
            "\nGot NMax1=$(NMax1) and NMax2=$(NMax2).",
        ),
    )
end

function Base.:+(
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
        return T3(pn1.coeffs[i] + pn2.coeffs[i])
    else
        return T3(pn2.coeffs[i])
    end
end

function Base.:*(
    pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}
) where {N1,N2,T1,T2,NMax1,NMax2}
    throw(
        ArgumentError(
            "`PNExpansion` multiplication is only defined for objects of the same PN order." *
            "\nGot NMax1=$(NMax1) and NMax2=$(NMax2).",
        ),
    )
end

function Base.:*(
    pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    if N1 > N2
        return pn2 * pn1
    else
        N3 = min(N1 + N2 - 1, NMax)
        PNExpansion(ntuple(i -> product_term(i, pn1, pn2), Val(N3)), NMax)
    end
end

function product_term(
    i, pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}
) where {N1,N2,T1,T2,NMax}
    T3 = promote_type(T1, T2)
    return sum(
        pn1.coeffs[j] * pn2.coeffs[i - j + 1] for
        j ∈ max(1, i - N2 + 1):min(i, N1), init ∈ zero(T3)
    )
end

Base.:/(p::PNExpansion, x::Number) = p * (1 / x)

function FastDifferentiation.derivative(
    pn_expansion::PNExpansion{N,T,NMax}, fd_node::FastDifferentiation.Node
) where {N,T,NMax}
    return PNExpansion(
        ntuple(i -> FastDifferentiation.derivative(pn_expansion[i], fd_node), Val(N)), NMax
    )
end

Base.Tuple(pn::PNExpansion) = pn.coeffs
SVector(pn::PNExpansion) = SVector(pn.coeffs)

function Base.:+(
    x::T1, term::PNTerm{T2,PNOrder,c⁻¹Exponent}
) where {T1<:Number,T2,PNOrder,c⁻¹Exponent}
    if c⁻¹exp(term) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " *
                "c⁻¹exp(term)=$(c⁻¹exp(term))." *
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    T = promote_type(T1, T2)
    N₀ = c⁻¹exp(term) + 1
    NMax = Int(2PNOrder + 1)
    N = min(N₀, NMax)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds coeffs[1] = x
    @inbounds if N₀ ≤ NMax
        coeffs[N₀] += term.coeff
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end
Base.:+(term::PNTerm, x::Number) = x + term

function Base.:-(term::PNTerm{T,PNOrder,c⁻¹Exponent}) where {T,PNOrder,c⁻¹Exponent}
    return PNTerm{T,PNOrder,c⁻¹Exponent}(-term.coeff)
end

function Base.:+(
    term1::PNTerm{T1,PNOrder,c⁻¹E1}, term2::PNTerm{T2,PNOrder,c⁻¹E2}
) where {T1,T2,PNOrder,c⁻¹E1,c⁻¹E2}
    if c⁻¹exp(term1) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " *
                "c⁻¹exp(term1)=$(c⁻¹exp(term1))." *
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    if c⁻¹exp(term2) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " *
                "c⁻¹exp(term2)=$(c⁻¹exp(term2))." *
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    T = promote_type(T1, T2)
    N1₀ = c⁻¹exp(term1) + 1
    N2₀ = c⁻¹exp(term2) + 1
    NMax = Int(2PNOrder + 1)
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

Base.:-(term1::PNTerm, term2::PNTerm) = term1 + (-term2)

function Base.:+(
    term::PNTerm{T1,PNOrder,c⁻¹E1}, expansion::PNExpansion{N2,T2,NMax2}
) where {T1,PNOrder,c⁻¹E1,N2,T2,NMax2}
    if c⁻¹exp(term) < 0
        throw(
            ArgumentError(
                "Cannot add a `PNTerm` with a negative exponent: " *
                "c⁻¹exp(term)=$(c⁻¹exp(term))." *
                "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
            ),
        )
    end
    N1 = c⁻¹exp(term) + 1
    NMax1 = Int(2PNOrder + 1)
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
Base.:+(expansion::PNExpansion, term::PNTerm) = term + expansion

Base.:-(term::PNTerm, x::Number) = term + (-x)
Base.:-(x::Number, term::PNTerm) = x + (-term)
Base.:-(term::PNTerm, expansion::PNExpansion) = term + (-expansion)
Base.:-(expansion::PNExpansion, term::PNTerm) = expansion + (-term)
Base.:-(x::Number, expansion::PNExpansion) = x + (-expansion)
Base.:-(expansion::PNExpansion, x::Number) = expansion + (-x)

function Base.:*(
    expansion::PNExpansion{N1,T1,NMax1}, term::PNTerm{T2,PNOrder,c⁻¹E2}
) where {N1,T1,NMax1,T2,PNOrder,c⁻¹E2}
    ΔN = c⁻¹exp(term)  # Note that ΔN may be negative!
    NMax2 = Int(2PNOrder + 1)
    NMax = min(NMax1, NMax2)
    N = min(max(N1, N1 + ΔN), NMax)

    # Check that no terms from expansion will be lost to negative PN orders
    @inbounds for i ∈ 1:min(max(0, -ΔN), N1)
        if !iszero(expansion[i])
            throw(
                ArgumentError(
                    "Cannot multiply `PNExpansion` by `PNTerm` with negative exponent: " *
                    "c⁻¹exp(term)=$(c⁻¹exp(term))." *
                    "\nResult will be a `PNExpansion`, which cannot store positive exponents.",
                ),
            )
        end
    end

    T = promote_type(T1, T2)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds for i ∈ max(1, 1 - ΔN):min(N1, N - ΔN)
        coeffs[i + ΔN] = expansion[i] * term.coeff
    end
    return PNExpansion{N,T,NMax}(Tuple(coeffs))
end
Base.:*(term::PNTerm, expansion::PNExpansion) = expansion * term
# (a, b, c, d, e, f, g) * (c⁻¹^2) = (0, 0, a, b, c, d, e)

Base.:/(expansion::PNExpansion, term::PNTerm) = expansion * inv(term)

@testitem "PNExpansion algebra" begin
    using Symbolics: @variables, simplify, substitute
    using PostNewtonian: PNExpansion

    for N1 ∈ 1:9
        for N2 ∈ 1:9
            for NMax ∈ max(N1, N2):(N1 + N2 + 3)
                @variables c⁻¹ x[1:N1] y[1:N2] z
                poly(e::PNExpansion) = sum(e[i] * c⁻¹^(i - 1) for i ∈ 1:length(e))
                eˣ = PNExpansion(tuple(x...), NMax)
                eʸ = PNExpansion(tuple(y...), NMax)

                # Test sums
                polysum = simplify(poly(eˣ + eʸ); expand=true)
                sumpoly = simplify(poly(eˣ) + poly(eʸ); expand=true)
                Δ = simplify(polysum - sumpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ + PNExpansion(tuple(z, x...), NMax + 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax + 1) + eˣ

                # Test products
                polyprod = simplify(poly(eˣ * eʸ); expand=true)
                prodpoly = simplify(
                    substitute(
                        simplify(poly(eˣ) * poly(eʸ); expand=true),
                        Dict([c⁻¹^n => 0 for n ∈ NMax:(2NMax + 3)]),
                    );
                    expand=true,
                )
                Δ = simplify(polyprod - prodpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ * PNExpansion(tuple(z, x...), NMax + 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax + 1) * eˣ
            end
        end
    end
end

@doc raw"""
    @pn_expansion [pnsystem] expansion

Mark `expansion` as a post-Newtonian expansion in powers of ``1/c``, to be truncated at the
order given by the `pnsystem`.

Note that we generally want this PN order to be interpreted as a *relative* PN order.  For
example, the expression for [gravitational-wave energy flux](@ref 𝓕) looks like
```math
\mathcal{F} = \frac{32c^5}{5G} \nu^2 \left(\frac{v}{c}\right)^{10}
\left[1 - \left(\frac{ν}{12} + \frac{3}{4}\right) \left(\frac{v}{c}\right)^2 + \ldots\right]
```
Here, the part in square brackets is the expansion in terms of *relative* PN order, with the
1 representing the 0-pN term, the ``(v/c)^2`` term being the 1-pN term, etc.  So in the
code, we apply `@pn_expansion` to the part in square brackets:
```julia
32c^5 / 5G * ν^2 * (v / c)^10 * @pn_expansion(1 - (ν/12 + 3/4) * (v/c)^2 + ...)
```
Note that the entire expression is written in a function that is modified by the
[`@pn_expression`](@ref) macro, which automatically inserts the `pnsystem` argument to this
`@pn_expansion` call, so you generally should not have to worry about it.

This macro gathers terms in `expansion` by the powers of ``1/c`` involved, zeroing out any
terms with powers of ``1/c`` higher than (twice) the `pnsystem`'s `PNOrder` parameter, and
combine the terms using the `PNExpansionReducer` specified in argument of the function that
includes this macro call.

The expansion and truncation are achieved by redefining `c` to be a
[`PNExpansionParameter(pnsystem)`](@ref PNExpansionParameter), which is just a
[`PNTerm`](@ref) with a coefficient of 1 and a `c⁻¹Exponent` of -1.  This redefinition
happens inside a `let` block created by this macro so that it doesn't interfere with any
`c`s on the outside.  For example, in the flux expression above, the `c`s in `32c^5 / 5G *
ν^2 * (v / c)^10` are outside the scope of the macro, even if we set `c=1` above that line,
the expansion will still be correct.
"""
@public macro pn_expansion(pnsystem, expr)
    return esc(quote
        let c = PNExpansionParameter($pnsystem)
            PNExpansionReducer($expr)
        end
    end)
end

@testitem "@pn_expansion" begin
    using PostNewtonian: @pn_expansion
    using MacroTools

    input = @macroexpand @pn_expansion pnsystem (1 - (ν/12 + 3/4) * (v/c)^2)
    output = quote
        let c = PNExpansionParameter(pnsystem)
            PNExpansionReducer(1 - (ν/12 + 3/4) * (v/c)^2)
        end
    end
    @test MacroTools.striplines(input) == MacroTools.striplines(output)
end
