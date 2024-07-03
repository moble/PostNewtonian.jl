# TODO: Remove N from PNExpansion; it should just always be MaxN

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
struct PNExpansion{N,T,MaxN}
    coeffs::NTuple{N,T}
end
#PNExpansion{N,T,MaxN}(coeffs::NTuple{N,T}) where {N,T,MaxN} = new{N,T,MaxN}(coeffs)

function PNExpansion(coeffs::NTuple{N,T}, ::Val{MaxN}) where {N,T,MaxN}
    # if MaxN < 1
    #     throw(ArgumentError("`MaxN` must be >0."))
    # end
    # if N > MaxN
    #     throw(ArgumentError("Length of `coeffs` is $N but must be <$MaxN."))
    # end
    PNExpansion{N,T,MaxN}(coeffs)
end

function PNExpansion(coeffs::NTuple{N,T}, MaxN::Int=typemax(Int) - 1) where {N,T}
    # if MaxN < 1
    #     throw(ArgumentError("`MaxN` must be >0."))
    # end
    # if N > MaxN
    #     throw(ArgumentError("Length of `coeffs` is $N but must be <$MaxN."))
    # end
    PNExpansion{N,T,MaxN}(coeffs)
end

pn_order(::PNExpansion{N,T,MaxN}) where {N,T,MaxN} = (MaxN - 1) // 2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = length(pn.coeffs)
Base.eltype(pn::PNExpansion) = eltype(pn.coeffs)

function Base.sum(pn_expansion::PNExpansion{N,T}) where {N,T}
    sum(pn_expansion[i] for i ∈ 1:N, init = zero(T))
end


function Base.:+(pn::PNExpansion{N,T1,MaxN}, x::T2) where {N,T1,MaxN,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(ntuple(i -> i == 1 ? pn[1] + x : T3(pn[i]), Val(N)), Val(MaxN))
end
Base.:+(x::T, pn::PNExpansion) where {T<:Number} = pn + x

function Base.:*(pn::PNExpansion{N,T1,MaxN}, x::T2) where {N,T1,MaxN,T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion((@. T3(pn.coeffs * x)), Val(MaxN))
end
Base.:*(x::T, pn::PNExpansion) where {T<:Number} = pn * x

function Base.:+(
    pn1::PNExpansion{N1,T1,MaxN1}, pn2::PNExpansion{N2,T2,MaxN2}
) where {N1,N2,T1,T2,MaxN1,MaxN2}
    error(
        "`PNExpansion` addition is only defined for objects of the same PN order."
        *
        "\nGot MaxN1=$(MaxN1) and MaxN2=$(MaxN2)."
    )
end

function Base.:+(
    pn1::PNExpansion{N1,T1,MaxN}, pn2::PNExpansion{N2,T2,MaxN}
) where {N1,N2,T1,T2,MaxN}
    if N1 > N2
        return pn2 + pn1
    else
        PNExpansion(ntuple(i -> sum_term(i, pn1, pn2), Val(N2)), Val(MaxN))
    end
end

function sum_term(
    i, pn1::PNExpansion{N1,T1,MaxN}, pn2::PNExpansion{N2,T2,MaxN}
) where {N1,N2,T1,T2,MaxN}
    T3 = promote_type(T1, T2)
    if i ≤ N1
        return T3(pn1.coeffs[i] + pn2.coeffs[i])
    else
        return T3(pn2.coeffs[i])
    end
end

function Base.:*(
    pn1::PNExpansion{N1,T1,MaxN1}, pn2::PNExpansion{N2,T2,MaxN2}
) where {N1,N2,T1,T2,MaxN1,MaxN2}
    error(
        "`PNExpansion` multiplication is only defined for objects of the same PN order."
        *
        "\nGot MaxN1=$(MaxN1) and MaxN2=$(MaxN2)."
    )
end

function Base.:*(
    pn1::PNExpansion{N1,T1,MaxN}, pn2::PNExpansion{N2,T2,MaxN}
) where {N1,N2,T1,T2,MaxN}
    if N1 > N2
        return pn2 * pn1
    else
        N3 = min(N1 + N2 - 1, MaxN)
        PNExpansion(ntuple(i -> product_term(i, pn1, pn2), Val(N3)), Val(MaxN))
    end
end

function product_term(
    i, pn1::PNExpansion{N1,T1,MaxN}, pn2::PNExpansion{N2,T2,MaxN}
) where {N1,N2,T1,T2,MaxN}
    T3 = promote_type(T1, T2)
    sum(
        pn1.coeffs[j] * pn2.coeffs[i-j+1] for j ∈ max(1, i - N2 + 1):min(i, N1),
        init = zero(T3)
    )
end

Base.:/(p::PNExpansion, x::Number) = p * (1 / x)



# TODO: Change PNExpansionParameter to PNTerm, which
#       - contains MaxExp in the type
#       - tracks the PN order of the term
#       - contains the value (coefficient) of this term
#       - can be exponentiated
#       - can divide a number
#       - can be multiplied by a number
#       - can be added to a number or another PNTerm, resulting in a PNExpansion

"""
    struct PNTerm{PNOrder}

This object represents a single term in a PNExpansion.  It has two fields: `exp`, which
represents the exponent of the PN expansion parameter ``1/c`` and thus will always be an
integer, usually positive; and `coeff`, which is the coefficient of the term.  The type
parameter `PNOrder` is a half-integer (just as in [`PNSystem`](@ref)s) representing the
order of the PN expansion.  If `exp > 2PNOrder`, then `coeff` will be zero.

"""
struct PNTerm{T,PNOrder}
    exp::Int
    coeff::T

    function PNTerm{T,PNOrder}(exp, coeff) where {T,PNOrder}
        if exp > 2PNOrder
            coeff = zero(coeff)
        end
        new{T,PNOrder}(exp, coeff)
    end
end
PNTerm(pn_order; exp=-1, coeff::T=1) where {T} = PNTerm{T,pn_order}(exp, coeff)
PNTerm(exp::Int, coeff::T, pn_order::Rational{Int}) where {T} = PNTerm{T,pn_order}(exp, coeff)

Base.length(pn::PNTerm) = 1
Base.eltype(pn::PNTerm{T}) where {T} = T

function Base.:inv(term::PNTerm{T,PNOrder}) where {T,PNOrder}
    PNTerm{T,PNOrder}(-term.exp, inv(term.coeff))
end

function Base.:^(term::PNTerm{T,PNOrder}, n::Int) where {T,PNOrder}
    coeff = term.coeff^n
    PNTerm{typeof(coeff),PNOrder}(term.exp * n, coeff)
end

function Base.:*(x::Number, term::PNTerm{T,PNOrder}) where {T,PNOrder}
    coeff = x * term.coeff
    PNTerm{typeof(coeff),PNOrder}(term.exp, coeff)
end
Base.:*(term::PNTerm, x::Number) = x * term

function Base.:/(term::PNTerm{T,PNOrder}, x::Number) where {T,PNOrder}
    coeff = term.coeff / x
    PNTerm{typeof(coeff),PNOrder}(term.exp, coeff)
end

function Base.:/(x::Number, term::PNTerm{T,PNOrder}) where {T,PNOrder}
    coeff = x / term.coeff
    PNTerm{typeof(coeff),PNOrder}(-term.exp, coeff)
end

function Base.:*(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    exp = term1.exp + term2.exp
    coeff = term1.coeff * term2.coeff
    PNTerm{typeof(coeff),PNOrder}(exp, coeff)
end

function Base.:/(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    exp = term1.exp - term2.exp
    coeff = term1.coeff / term2.coeff
    PNTerm{typeof(coeff),PNOrder}(exp, coeff)
end

function Base.:+(x::T1, term::PNTerm{T2,PNOrder}) where {T1<:Number,T2,PNOrder}
    if term.exp < 0
        throw(ArgumentError(
            "Cannot *add* a `PNTerm` with a negative exponent: "
            * "term.exp=$(term.exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    MaxN = Int(2PNOrder + 1)
    # coeffs = ntuple(
    #     i -> (
    #         (i==1 ? T(x) : zero(T))
    #         +
    #         (term.exp == i - 1 ? T(term.coeff) : zero(T))
    #     ),
    #     Val(MaxN)
    # )
    # PNExpansion{MaxN,T,MaxN}(coeffs)
    coeffs = MVector{MaxN, T}(undef)
    coeffs .= zero(T)
    @inbounds coeffs[1] = x
    @inbounds if term.exp < MaxN
        coeffs[term.exp + 1] += term.coeff
    end
    PNExpansion{MaxN,T,MaxN}(Tuple(coeffs))
end
Base.:+(term::PNTerm{T1,PNOrder}, x::T2) where {T1,T2<:Number,PNOrder} = x + term

function Base.:+(term1::PNTerm{T1,PNOrder}, term2::PNTerm{T2,PNOrder}) where {T1,T2,PNOrder}
    if term1.exp < 0
        throw(ArgumentError(
            "Cannot *add* a `PNTerm` with a negative exponent: "
            * "term1.exp=$(term1.exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    if term2.exp < 0
        throw(ArgumentError(
            "Cannot *add* a `PNTerm` with a negative exponent: "
            * "term2.exp=$(term2.exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    MaxN = Int(2PNOrder + 1)
    coeffs = MVector{MaxN, T}(undef)
    coeffs .= zero(T)
    @inbounds if term1.exp < MaxN
        coeffs[term1.exp + 1] += term1.coeff
    end
    @inbounds if term2.exp < MaxN
        coeffs[term2.exp + 1] += term2.coeff
    end
    PNExpansion{MaxN,T,MaxN}(Tuple(coeffs))
    # coeffs = ntuple(
    #     i -> (
    #         (term1.exp == i + 1 ? T(term1.coeff) : zero(T))
    #         +
    #         (term2.exp == i + 1 ? T(term2.coeff) : zero(T))
    #     ),
    #     MaxN
    # )
    # PNExpansion{MaxN,T,MaxN}(coeffs)
end

function Base.:+(term::PNTerm{T1,PNOrder}, expansion::PNExpansion{N,T2,MaxN}) where {T1,T2,N,MaxN,PNOrder}
    if term.exp < 0
        throw(ArgumentError(
            "Cannot *add* a `PNTerm` with a negative exponent: "
            * "term.exp=$(term.exp)."
            * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
        ))
    end
    T = promote_type(T1, T2)
    coeffs = MVector{MaxN, T}(undef)
    coeffs .= zero(T)
    @inbounds if term.exp < MaxN
        coeffs[term.exp + 1] += term.coeff
    end
    @inbounds for i ∈ 1:N
        if i ≤ MaxN
            coeffs[i] += expansion[i]
        end
    end
    PNExpansion{MaxN,T,MaxN}(Tuple(coeffs))
    # coeffs = ntuple(
    #     i -> (
    #         (term.exp == i + 1 ? T(term.coeff) : zero(T))
    #         +
    #         (i ≤ N ? T(expansion[i]) : zero(T))
    #     ),
    #     MaxN
    # )
    # PNExpansion{MaxN,T,MaxN}(coeffs)
end
Base.:+(expansion::PNExpansion, term::PNTerm) = term + expansion

function Base.:/(expansion::PNExpansion{N,T1,MaxN}, term::PNTerm{T2,PNOrder}) where {T1,T2,N,MaxN,PNOrder}
    T = promote_type(T1, T2)
    coeffs = MVector{MaxN, T}(undef)
    coeffs .= zero(T)
    @inbounds if term.exp < MaxN
        coeffs[term.exp + 1] += term.coeff
    end
    @inbounds for i ∈ 1:N
        if i ≤ MaxN
            coeffs[i] = expansion[i] / term.coeff
        end
    end
    PNExpansion{MaxN,T,MaxN}(Tuple(coeffs))
end



function PNExpansionParameter(::PNSystem{ST, PNOrder}) where {ST,PNOrder}
    PNTerm{eltype(ST), PNOrder}(-1, one(eltype(ST)))
end



# """
#     struct PNExpansionParameter

# A struct representing a Post-Newtonian expansion parameter. It has fields `exp` and `max_N`.

# ## Fields
# - `exp`: The exponent of the expansion parameter.
# - `max_N`: The maximum order of the expansion.

# ## Constructors
# - `PNExpansionParameter()`: Constructs a `PNExpansionParameter` with default values for `exp` and `max_N`.
# - `PNExpansionParameter(exp, max_N)`: Constructs a `PNExpansionParameter` with the given values for `exp` and `max_N`.
# - `PNExpansionParameter(; exp=1, max_N=typemax(Int)-1)`: Constructs a `PNExpansionParameter` with the given values for `exp` and `max_N`.

# ## Methods
# - `Base.:^(::PNExpansionParameter, n::Int)`: Raises a `PNExpansionParameter` to the power of an integer `n`.
# - `Base.:/(x, ::PNExpansionParameter)`: Divides a value `x` by a `PNExpansionParameter`.
# """
# struct PNExpansionParameter{MaxN}
#     exp::Int

#     PNExpansionParameter{MaxN}(exp=1) where {MaxN} = new{MaxN}(exp)
#     # function PNExpansionParameter(; exp=1, max_N=typemax(Int)-1, pn_order=(max_N-1)//2)
#     #     if pn_order != (max_N-1)//2
#     #         max_N = Int(2pn_order + 1)
#     #     end
#     #     if !isinteger(max_N) || max_N < 1
#     #         throw(ArgumentError("`max_N` must be a positive integer; got $max_N."))
#     #     end
#     #     new(exp, max_N)
#     # end
# end

# function Base.:^(c::PNExpansionParameter{MaxN}, n::Int) where {MaxN}
#     PNExpansionParameter{MaxN}(n * c.exp)
# end

# function Base.:/(x, c::PNExpansionParameter{MaxN}) where {MaxN}
#     if c.exp < 0
#         throw(ArgumentError(
#             "Cannot *divide* by a `PNExpansionParameter` raised to a negative power: "
#             * "c.exp=$(c.exp)."
#             * "\nResult will be a `PNExpansion`, which cannot store positive exponents."
#         ))
#     end
#     if 1 + c.exp > MaxN
#         PNExpansion(NTuple{1,typeof(x)}(zero(x)), MaxN)
#     else
#         PNExpansion(ntuple(i -> i == 1 + c.exp ? x : zero(x), 1 + c.exp), MaxN)
#     end
# end

# const c = PNExpansionParameter{typemax(Int) - 1}()#1, typemax(Int) - 1)
