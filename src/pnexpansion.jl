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
struct PNExpansion{N, T, MaxN}
    coeffs::NTuple{N, T}

    function PNExpansion(coeffs::NTuple{N, T}, MaxN::Int=typemax(Int)-1) where {N, T}
        if MaxN < 1
            throw(ArgumentError("`MaxN` must be >0."))
        end
        if N > MaxN
            throw(ArgumentError("Length of `coeffs` is $N but must be <$MaxN."))
        end
        new{N, T, MaxN}(coeffs)
    end
end
pn_order(::PNExpansion{N, T, MaxN}) where {N, T, MaxN} = (MaxN-1) // 2

Base.getindex(pn::PNExpansion, i::Int) = pn.coeffs[i]
Base.length(pn::PNExpansion) = length(pn.coeffs)
Base.show(io::IO, pn::PNExpansion) = print(io, "PNExpansion($(pn.coeffs))")

function Base.:+(pn::PNExpansion{N, T1, MaxN}, x::T2) where {N, T1, MaxN, T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(ntuple(i-> i==1 ? pn[1]+x : T3(pn[i]), Val(N)), MaxN)
end
Base.:+(x::T, pn::PNExpansion) where {T<:Number} = pn + x

function Base.:*(pn::PNExpansion{N, T1, MaxN}, x::T2) where {N, T1, MaxN, T2<:Number}
    T3 = promote_type(T1, T2)
    PNExpansion(@. T3(pn * x), MaxN)
end
Base.:*(x::T, pn::PNExpansion) where {T<:Number} = pn * x

function Base.:+(
    pn1::PNExpansion{N1, T1, MaxN1}, pn2::PNExpansion{N2, T2, MaxN2}
) where {N1, N2, T1, T2, MaxN1, MaxN2}
    error(
        "`PNExpansion` addition is only defined for objects of the same PN order."
        * "\nGot MaxN1=$(MaxN1) and MaxN2=$(MaxN2)."
    )
end

function Base.:+(
    pn1::PNExpansion{N1, T1, MaxN}, pn2::PNExpansion{N2, T2, MaxN}
) where {N1, N2, T1, T2, MaxN}
    if N1 > N2
        return pn2 + pn1
    else
        PNExpansion(ntuple(i->sum_term(i, pn1, pn2), Val(N2)), MaxN)
    end
end

function sum_term(
    i, pn1::PNExpansion{N1, T1, MaxN}, pn2::PNExpansion{N2, T2, MaxN}
) where {N1, N2, T1, T2, MaxN}
    T3 = promote_type(T1, T2)
    if i ≤ N1
        return T3(pn1.coeffs[i] + pn2.coeffs[i])
    else
        return T3(pn2.coeffs[i])
    end
end

function Base.:*(
    pn1::PNExpansion{N1, T1, MaxN1}, pn2::PNExpansion{N2, T2, MaxN2}
) where {N1, N2, T1, T2, MaxN1, MaxN2}
    error(
        "`PNExpansion` multiplication is only defined for objects of the same PN order."
        * "\nGot MaxN1=$(MaxN1) and MaxN2=$(MaxN2)."
    )
end

function Base.:*(
    pn1::PNExpansion{N1, T1, MaxN}, pn2::PNExpansion{N2, T2, MaxN}
) where {N1, N2, T1, T2, MaxN}
    if N1 > N2
        return pn2 * pn1
    else
        N3 = min(N1+N2-1, MaxN)
        PNExpansion(ntuple(i->product_term(i, pn1, pn2), Val(N3)), MaxN)
    end
end

function product_term(
    i, pn1::PNExpansion{N1, T1, MaxN}, pn2::PNExpansion{N2, T2, MaxN}
) where {N1, N2, T1, T2, MaxN}
    T3 = promote_type(T1, T2)
    sum(
        pn1.coeffs[j]*pn2.coeffs[i-j+1] for j ∈ max(1,i-N2+1):min(i,N1),
        init=zero(T3)
    )
end

Base.:/(p::PNExpansion, x::Number) = p * (1/x)



"""
    struct PNExpansionParameter{Exp, MaxN}

A struct representing a Post-Newtonian expansion parameter. It is parameterized by `Exp` and `MaxN`.

## Fields
- `Exp`: The exponent of the expansion parameter.
- `MaxN`: The maximum order of the expansion.

## Constructors
- `PNExpansionParameter{Exp, MaxN}() where {Exp, MaxN}`: Constructs a `PNExpansionParameter` with default values for `Exp` and `MaxN`.
- `PNExpansionParameter(exp, maxN)`: Constructs a `PNExpansionParameter` with the given values for `Exp` and `MaxN`.
- `PNExpansionParameter(; exp=1, pnorder=(typemax(Int)-2)//2)`: Constructs a `PNExpansionParameter` with the given values for `exp` and `pnorder`.

## Methods
- `Base.:^(::PNExpansionParameter{Exp, MaxN}, n::Int) where {Exp, MaxN}`: Raises a `PNExpansionParameter` to the power of an integer `n`.
- `Base.:/(x, ::PNExpansionParameter{Exp, MaxN}) where {Exp, MaxN}`: Divides a value `x` by a `PNExpansionParameter`.
"""
struct PNExpansionParameter{Exp, MaxN}
    PNExpansionParameter{Exp, MaxN}() where {Exp, MaxN} = new{Exp, MaxN}()
    PNExpansionParameter(exp, maxN) = new{exp, maxN}()
    function PNExpansionParameter(; exp=1, pnorder=(typemax(Int)-2)//2)
        MaxN = 2pnorder + 1
        if !isinteger(MaxN)
            throw(ArgumentError("`pnorder` must be a half-integer; got $pnorder."))
        end
        new{exp, Int(MaxN)}()
    end
end

function Base.:^(::PNExpansionParameter{Exp, MaxN}, n::Int) where {Exp, MaxN}
    PNExpansionParameter(Exp*n, MaxN)
    # PNExpansionParameter{Exp*n, MaxN}()
end

function Base.:/(x, ::PNExpansionParameter{Exp, MaxN}) where {Exp, MaxN}
    if Exp < 0
        throw(ArgumentError(
            "Cannot divide by a PNExpansionParameter raised to a negative power."
        ))
    end
    if 1+Exp > MaxN
        PNExpansion(NTuple{1, typeof(x)}(zero(x)), MaxN)
    else
        PNExpansion(ntuple(i->i==1+Exp ? x : zero(x), Val(1+Exp)), MaxN)
    end
end

const c = PNExpansionParameter{1, typemax(Int)-1}()
