# This a utility that allow us to interoperate with FastDifferentiation.Node and other
# Number types.
function _efficient_vector(::Val{N}, ::Val{T}) where {N,T}
    if isbitstype(T)
        MVector{N, T}(undef)
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
tuple of coefficients, and `T` is the type of the coefficients.  The `NMax` parameter is the
maximum number of elements, related to the usual PN order by
```math
\text{pn_order} = \frac{\texttt{NMax}-1} {2}.
```
The `N` parameter is not related to the PN order; it is just used by Julia to know how many
elements are currently in the coefficients, but is required to be 1 ≤ N ≤ NMax.

"""
struct PNExpansion{N,T,NMax}
    coeffs::NTuple{N,T}

    function PNExpansion{N,T,NMax}(coeffs) where {N,T,NMax}
        if N < 1
            throw(ArgumentError("`N=$N` must be >0."))
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be <`NMax=$NMax`."))
        end
        new{N,T,NMax}(coeffs)
    end
    function PNExpansion(coeffs::NTuple{N,T}, NMax) where {N,T}
        if N < 1
            throw(ArgumentError("`N=$N` must be >0."))
        end
        if N > NMax
            throw(ArgumentError("`N=$N` must be <`NMax=$NMax`."))
        end
        new{N,T,NMax}(coeffs)
    end
end

pn_order(::PNExpansion{N,T,NMax}) where {N,T,NMax} = (NMax - 1) // 2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = length(pn.coeffs)
Base.eltype(pn::PNExpansion) = eltype(pn.coeffs)

function Base.sum(pn_expansion::PNExpansion{N,T,NMax}) where {N,T,NMax}
    sum(pn_expansion[i] for i ∈ 1:N, init = zero(T))
end

function Base.:+(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(ntuple(i -> i == 1 ? pn[1] + x : T3(pn[i]), Val(N)), NMax)
end
Base.:+(x::T, pn::PNExpansion) where {T<:Number} = pn + x

function Base.:-(pn::PNExpansion{N,T,NMax}) where {N,T,NMax}
    PNExpansion{N,T,NMax}((-).(pn.coeffs))
end

function Base.:*(pn::PNExpansion{N,T1,NMax}, x::T2) where {N,T1,NMax,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion{N,T3,NMax}(@. T3(pn.coeffs * x))
end
Base.:*(x::T, pn::PNExpansion) where {T<:Number} = pn * x

function Base.:+(pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}) where
{N1,N2,T1,T2,NMax1,NMax2}
    throw(ArgumentError(
        "`PNExpansion` addition is only defined for objects of the same PN order."
        *
        "\nGot NMax1=$(NMax1) and NMax2=$(NMax2)."
    ))
end

function Base.:+(pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}) where
{N1,N2,T1,T2,NMax}
    if N1 > N2
        return pn2 + pn1
    else
        PNExpansion(ntuple(i -> sum_term(i, pn1, pn2), Val(N2)), NMax)
    end
end

function sum_term(i, pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}) where
{N1,N2,T1,T2,NMax}
    T3 = promote_type(T1, T2)
    if i ≤ N1
        return T3(pn1.coeffs[i] + pn2.coeffs[i])
    else
        return T3(pn2.coeffs[i])
    end
end

function Base.:*(pn1::PNExpansion{N1,T1,NMax1}, pn2::PNExpansion{N2,T2,NMax2}) where
{N1,N2,T1,T2,NMax1,NMax2}
    throw(ArgumentError(
        "`PNExpansion` multiplication is only defined for objects of the same PN order."
        *
        "\nGot NMax1=$(NMax1) and NMax2=$(NMax2)."
    ))
end

function Base.:*(pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}) where
{N1,N2,T1,T2,NMax}
    if N1 > N2
        return pn2 * pn1
    else
        N3 = min(N1 + N2 - 1, NMax)
        PNExpansion(ntuple(i -> product_term(i, pn1, pn2), Val(N3)), NMax)
    end
end

function product_term(i, pn1::PNExpansion{N1,T1,NMax}, pn2::PNExpansion{N2,T2,NMax}) where
{N1,N2,T1,T2,NMax}
    T3 = promote_type(T1, T2)
    sum(
        pn1.coeffs[j] * pn2.coeffs[i-j+1] for j ∈ max(1, i - N2 + 1):min(i, N1),
        init = zero(T3)
    )
end

Base.:/(p::PNExpansion, x::Number) = p * (1 / x)

function FastDifferentiation.derivative(
    pn_expansion::PNExpansion{N,T,NMax},
    fd_node::FastDifferentiation.Node
) where {N,T,NMax}
    PNExpansion(
        ntuple(
            i -> FastDifferentiation.derivative(pn_expansion[i], fd_node),
            Val(N)
        ),
        NMax
    )
end

Base.Tuple(pn::PNExpansion) = pn.coeffs
SVector(pn::PNExpansion) = SVector(pn.coeffs)


"""
    struct PNTerm{T,PNOrder}

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

function Base.sum(pn::PNTerm)
    pn.coeff
end

function Base.:+(pn::PNTerm)
    pn
end

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
    N₀ = term.c⁻¹exp + 1
    NMax = Int(2PNOrder + 1)
    N = min(N₀, NMax)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds coeffs[1] = x
    @inbounds if N₀ ≤ NMax
        coeffs[N₀] += term.coeff
    end
    PNExpansion{N,T,NMax}(Tuple(coeffs))
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
    N1₀ = term1.c⁻¹exp + 1
    N2₀ = term2.c⁻¹exp + 1
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
    PNExpansion{N,T,NMax}(Tuple(coeffs))
end
Base.:-(term1::PNTerm, term2::PNTerm) = term1 + (-term2)

function Base.:+(term::PNTerm{T1,PNOrder}, expansion::PNExpansion{N2,T2,NMax2}) where
{T1,PNOrder,N2,T2,NMax2}
    if term.c⁻¹exp < 0
        throw(ArgumentError(
            "Cannot add a `PNTerm` with a negative exponent: "
            * "term.c⁻¹exp=$(term.c⁻¹exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    N1 = term.c⁻¹exp + 1
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
    PNExpansion{N,T,NMax}(Tuple(coeffs))
end
Base.:+(expansion::PNExpansion, term::PNTerm) = term + expansion

Base.:-(term::PNTerm, x::Number) = term + (-x)
Base.:-(x::Number, term::PNTerm) = x + (-term)
Base.:-(term::PNTerm, expansion::PNExpansion) = term + (-expansion)
Base.:-(expansion::PNExpansion, term::PNTerm) = expansion + (-term)
Base.:-(x::Number, expansion::PNExpansion) = x + (-expansion)
Base.:-(expansion::PNExpansion, x::Number) = expansion + (-x)

function Base.:*(expansion::PNExpansion{N1,T1,NMax1}, term::PNTerm{T2,PNOrder}) where
{N1,T1,NMax1,T2,PNOrder}
    ΔN = term.c⁻¹exp  # Note that ΔN may be negative!
    NMax2 = Int(2PNOrder + 1)
    NMax = min(NMax1, NMax2)
    N = min(max(N1, N1+ΔN), NMax)

    # Check that no terms from expansion will be lost to negative PN orders
    @inbounds for i ∈ 1:min(max(0,-ΔN), N1)
        if !iszero(expansion[i])
            throw(ArgumentError(
                "Cannot multiply `PNExpansion` by `PNTerm` with negative exponent: "
                * "term.c⁻¹exp=$(term.c⁻¹exp)."
                * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
            ))
        end
    end

    T = promote_type(T1, T2)
    coeffs = _efficient_vector(Val(N), Val(T))
    coeffs .= zero(T)
    @inbounds for i ∈ max(1,1-ΔN):min(N1,N-ΔN)
        coeffs[i+ΔN] = expansion[i] * term.coeff
    end
    PNExpansion{N,T,NMax}(Tuple(coeffs))
end
Base.:*(term::PNTerm, expansion::PNExpansion) = expansion * term
# (a, b, c, d, e, f, g) * (c⁻¹^2) = (0, 0, a, b, c, d, e)

Base.:/(expansion::PNExpansion, term::PNTerm) = expansion * inv(term)


function PNExpansionParameter(::PNSystem{ST, PNOrder}) where {ST,PNOrder}
    PNTerm{eltype(ST), PNOrder}(-1, one(eltype(ST)))
end
