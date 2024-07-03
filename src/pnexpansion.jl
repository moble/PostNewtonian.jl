@doc raw"""

This object can be multiplied by a scalar or another `PNExpansion` object, and contains a
tuple of coefficients.  The coefficients are stored in the order of the expansion, with the
zeroth-order coefficient first.  Multiplication by a scalar just multiplies each of the
elements.  Multiplication by another `PNExpansion` object is done by convolving the two
tuples of coefficients.

Blanchet (2014) defines the post-Newtonian expansion parameter as follows:

> This parameter represents essentially a slow motion estimate ``𝜖 ∼ 𝑣/𝑐``, where ``𝑣``
> denotes a typical internal velocity.  By a slight abuse of notation, following
> Chandrasekhar et al. [...], we shall henceforth write formally ``𝜖 ≡ 1/𝑐``, even though
> ``𝜖`` is dimensionless whereas ``𝑐`` has the dimension of a velocity. Thus, ``1/𝑐 ≪ 1``
> in the case of post-Newtonian sources. The small post-Newtonian remainders will be denoted
> ``𝒪(1/𝑐^𝑛)``. Furthermore, [...] we shall refer to a small post-Newtonian term with
> formal order ``𝒪(1/𝑐^𝑛)`` relative to the Newtonian acceleration in the equations of
> motion, as ``\frac{𝑛}{2}\text{PN}``.

Therefore, we consider the coefficients of the `PNExpansion` to be a polynomial in ``1/𝑐``.
Here, the type parameter `N` corresponds to the number of elements actually present in the
tuple of coefficients, and `T` is the type of the coefficients.  The `MaxN` parameter is the
maximum order of the expansion, related to the usual PN order by
```math
\text{pn_order} = \frac{\texttt{MaxN}-1} {2}.
```
The `N` parameter is not related to the PN order; it is just used by Julia to know how many
elements are currently in the coefficients.

"""
struct PNExpansion{N,T}
    coeffs::NTuple{N,T}

    PNExpansion{N,T}(coeffs) where {N,T} = new{N,T}(coeffs)
    function PNExpansion(coeffs::NTuple{N,T}) where {N,T}
        if N < 1
            throw(ArgumentError("`N` must be >0."))
        end
        new{N,T}(coeffs)
    end
end

pn_order(::PNExpansion{N,T}) where {N,T} = (N - 1) // 2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = length(pn.coeffs)
Base.eltype(pn::PNExpansion) = eltype(pn.coeffs)

function Base.sum(pn_expansion::PNExpansion{N,T}) where {N,T}
    sum(pn_expansion[i] for i ∈ 1:N, init = zero(T))
end

function Base.:+(pn::PNExpansion{N,T1}, x::T2) where {N,T1,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(ntuple(i -> i == 1 ? pn[1] + x : T3(pn[i]), Val(N)))
end
Base.:+(x::T, pn::PNExpansion) where {T<:Number} = pn + x

function Base.:-(pn::PNExpansion)
    PNExpansion((-).(pn.coeffs))
end

function Base.:*(pn::PNExpansion{N,T1}, x::T2) where {N,T1,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(@. T3(pn.coeffs * x))
end
Base.:*(x::T, pn::PNExpansion) where {T<:Number} = pn * x

function Base.:+(pn1::PNExpansion{N1,T1}, pn2::PNExpansion{N2,T2}) where {N1,N2,T1,T2}
    error(
        "`PNExpansion` addition is only defined for objects of the same PN order."
        *
        "\nGot N1=$(N1) and N2=$(N2)."
    )
end

function Base.:+(pn1::PNExpansion{N,T1}, pn2::PNExpansion{N,T2}) where {N,T1,T2}
    PNExpansion(ntuple(i -> pn1.coeffs[i] + pn2.coeffs[i], Val(N)))
end

function Base.:*(pn1::PNExpansion{N1,T1}, pn2::PNExpansion{N2,T2}) where {N1,N2,T1,T2}
    error(
        "`PNExpansion` multiplication is only defined for objects of the same PN order."
        *
        "\nGot N1=$(N1) and N2=$(N2)."
    )
end

function Base.:*(pn1::PNExpansion{N,T1}, pn2::PNExpansion{N,T2}) where {N,T1,T2}
    PNExpansion(ntuple(i -> product_term(i, pn1, pn2), Val(N)))
end

function product_term(i, pn1::PNExpansion{N,T1}, pn2::PNExpansion{N,T2}) where {N,T1,T2}
    T3 = promote_type(T1, T2)
    sum(
        pn1.coeffs[j] * pn2.coeffs[i-j+1] for j ∈ max(1, i - N + 1):min(i, N),
        init = zero(T3)
    )
end

Base.:/(p::PNExpansion, x::Number) = p * (1 / x)


"""
    struct PNTerm{PNOrder}

This object represents a single term in a PNExpansion.  It has two fields: `c⁻¹exp`, which
represents the exponent of the PN expansion parameter ``1/c`` and thus will always be an
integer, usually positive; and `coeff`, which is the coefficient of the term.  The type
parameter `PNOrder` is a half-integer (just as in [`PNSystem`](@ref)s) representing the
order of the PN expansion.  If `c⁻¹exp > 2PNOrder`, then `coeff` will be zero.

"""
struct PNTerm{T,PNOrder}
    c⁻¹exp::Int
    coeff::T

    function PNTerm{T,PNOrder}(c⁻¹exp, coeff) where {T,PNOrder}
        if c⁻¹exp > 2PNOrder
            coeff = zero(coeff)
        end
        new{T,PNOrder}(c⁻¹exp, coeff)
    end
end

Base.length(pn::PNTerm) = 1
Base.eltype(pn::PNTerm{T}) where {T} = T

function Base.inv(term::PNTerm{T,PNOrder}) where {T,PNOrder}
    PNTerm{T,PNOrder}(-term.c⁻¹exp, inv(term.coeff))
end

function Base.:^(term::PNTerm{T,PNOrder}, n::Int) where {T,PNOrder}
    coeff = term.coeff^n
    PNTerm{typeof(coeff),PNOrder}(term.c⁻¹exp * n, coeff)
end

function Base.:*(x::Number, term::PNTerm{T,PNOrder}) where {T,PNOrder}
    coeff = x * term.coeff
    PNTerm{typeof(coeff),PNOrder}(term.c⁻¹exp, coeff)
end
Base.:*(term::PNTerm, x::Number) = x * term

function Base.:/(term::PNTerm{T,PNOrder}, x::Number) where {T,PNOrder}
    coeff = term.coeff / x
    PNTerm{typeof(coeff),PNOrder}(term.c⁻¹exp, coeff)
end

function Base.:/(x::Number, term::PNTerm{T,PNOrder}) where {T,PNOrder}
    coeff = x / term.coeff
    PNTerm{typeof(coeff),PNOrder}(-term.c⁻¹exp, coeff)
end

function Base.:+(x::T1, term::PNTerm{T2,PNOrder}) where {T1<:Number,T2,PNOrder}
    if term.c⁻¹exp < 0
        throw(ArgumentError(
            "Cannot add a `PNTerm` with a negative exponent: "
            * "term.c⁻¹exp=$(term.c⁻¹exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    N = Int(2PNOrder + 1)
    coeffs = MVector{N, T}(undef)
    coeffs .= zero(T)
    @inbounds coeffs[1] = x
    @inbounds if term.c⁻¹exp < N
        coeffs[term.c⁻¹exp + 1] += term.coeff
    end
    PNExpansion{N,T}(Tuple(coeffs))
end
Base.:+(term::PNTerm{T1,PNOrder}, x::T2) where {T1,T2<:Number,PNOrder} = x + term

function Base.:-(term::PNTerm{T,PNOrder}) where {T,PNOrder}
    PNTerm{T,PNOrder}(term.c⁻¹exp, -term.coeff)
end

function Base.:*(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    c⁻¹exp = term1.c⁻¹exp + term2.c⁻¹exp
    coeff = term1.coeff * term2.coeff
    PNTerm{typeof(coeff),PNOrder}(c⁻¹exp, coeff)
end

function Base.:/(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    c⁻¹exp = term1.c⁻¹exp - term2.c⁻¹exp
    coeff = term1.coeff / term2.coeff
    PNTerm{typeof(coeff),PNOrder}(c⁻¹exp, coeff)
end

function Base.:+(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    if term1.c⁻¹exp < 0
        throw(ArgumentError(
            "Cannot add a `PNTerm` with a negative exponent: "
            * "term1.c⁻¹exp=$(term1.c⁻¹exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    if term2.c⁻¹exp < 0
        throw(ArgumentError(
            "Cannot add a `PNTerm` with a negative exponent: "
            * "term2.c⁻¹exp=$(term2.c⁻¹exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    N = Int(2PNOrder + 1)
    coeffs = MVector{N, T}(undef)
    coeffs .= zero(T)
    @inbounds if term1.c⁻¹exp < N
        coeffs[term1.c⁻¹exp + 1] += term1.coeff
    end
    @inbounds if term2.c⁻¹exp < N
        coeffs[term2.c⁻¹exp + 1] += term2.coeff
    end
    PNExpansion{N,T}(Tuple(coeffs))
end
Base.:-(term1::PNTerm, term2::PNTerm) = term1 + (-term2)

function Base.:+(term::PNTerm{T1,PNOrder}, expansion::PNExpansion{N,T2}) where {T1,T2,N,PNOrder}
    if term.c⁻¹exp < 0
        throw(ArgumentError(
            "Cannot add a `PNTerm` with a negative exponent: "
            * "term.c⁻¹exp=$(term.c⁻¹exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    coeffs = MVector{N, T}(undef)
    coeffs .= zero(T)
    @inbounds if term.c⁻¹exp < N
        coeffs[term.c⁻¹exp + 1] += term.coeff
    end
    @inbounds for i ∈ 1:N
        if i ≤ N
            coeffs[i] += expansion[i]
        end
    end
    PNExpansion{N,T}(Tuple(coeffs))
end
Base.:+(expansion::PNExpansion, term::PNTerm) = term + expansion

Base.:-(term::PNTerm, x::Number) = term + (-x)
Base.:-(x::Number, term::PNTerm) = x + (-term)
Base.:-(term::PNTerm, expansion::PNExpansion) = term + (-expansion)
Base.:-(expansion::PNExpansion, term::PNTerm) = expansion + (-term)
Base.:-(x::Number, expansion::PNExpansion) = x + (-expansion)
Base.:-(expansion::PNExpansion, x::Number) = expansion + (-x)

function Base.:/(expansion::PNExpansion{N,T1}, term::PNTerm{T2,PNOrder}) where {N,T1,T2,PNOrder}
    T = promote_type(T1, T2)
    coeffs = MVector{N, T}(undef)
    coeffs .= zero(T)
    @inbounds for i ∈ max(1,1+term.c⁻¹exp):min(N,N+term.c⁻¹exp)
        coeffs[i - term.c⁻¹exp] = expansion[i] / term.coeff
    end
    PNExpansion{N,T}(Tuple(coeffs))
end


function PNExpansionParameter(::PNSystem{ST, PNOrder}) where {ST,PNOrder}
    PNTerm{eltype(ST), PNOrder}(-1, one(eltype(ST)))
end
