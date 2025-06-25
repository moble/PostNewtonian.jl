"""
    PNSystem{NT, PNOrder, ST}

Base type for all PN systems, such as `BBH`, `BHNS`, and `NSNS`.

These objects encode all essential properties of the binary, including its current state.
As such, they can be used as inputs to the various [fundamental](@ref Fundamental-variables)
and [derived variables](@ref Derived-variables), as well as [PN expressions](@ref) and
[dynamics](@ref Dynamics) functions.

The parameter `NT` is the number type of the system, such as `Float64` or `Dual{SomeTag,
Float64, 7}`.  The parameter `ST <: DenseVector{NT}` is the type returned by the `state`
function, which probably just returns the `state` vector stored in the concrete subtype.  As
such, this will probably be `MVector{N, NT}` or `SVector{N, NT}`, where `N` is the number of
elements in the state.  `PNOrder` is a `Rational` giving the order to which PN expansions
should be carried out when using the given object.
"""
abstract type PNSystem{NT,PNOrder,ST<:DenseVector{NT}} <: DenseVector{NT} end

"""
    state(pnsystem::PNSystem)

Return the state vector of `pnsystem`, which is a vector of fundamental variables for the
given PN system.

Note that the built-in `PNSystem` subtypes have a `state` field that is a vector, so this
function will just return that vector.  However, that may not always be true for
user-defined subtypes.
"""
function state(::T) where {T<:PNSystem}
    error("`state` is not yet defined for PNSystem subtype `$T`.")
end
Base.vec(pnsystem::PNSystem) = state(pnsystem)

Base.eltype(::Type{PNT}) where {NT,PNT<:PNSystem{NT}} = NT
Base.one(::Type{PNT}) where {PNT<:PNSystem} = one(eltype(PNT))
Base.one(x::T) where {T<:PNSystem} = one(T)
Base.zero(::Type{PNT}) where {PNT<:PNSystem} = zero(eltype(PNT))
Base.zero(x::T) where {T<:PNSystem} = zero(T)
Base.float(::Type{PNT}) where {PNT<:PNSystem} = float(eltype(PNT))
Base.float(x::T) where {T<:PNSystem} = float(T)

"""
    symbols(pnsystem::PNSystem)
    symbols(::Type{<:PNSystem})
    ascii_symbols(pnsystem::PNSystem)
    ascii_symbols(::Type{<:PNSystem})

Return a Tuple of symbols corresponding to the variables tracked by `pnsystem`, in the order
in which they are stored in the `state` vector.

The `ascii_symbols` function returns those symbols in ASCII form, enabling interaction with
external systems (e.g., Python) that do not support many Unicode symbols.

```jldoctest
julia> using PostNewtonian: BBH

julia> pnsystem = BBH(randn(14); PNOrder=7//2);

julia> symbols(pnsystem)
(:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)

julia> ascii_symbols(pnsystem)
(:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
```
"""
symbols(pnsystem::PNSystem) = symbols(typeof(pnsystem))
function symbols(::Type{T}) where {T<:PNSystem}
    error("`symbols` is not yet defined for PNSystem subtype `$T`.")
end
ascii_symbols(pnsystem::PNSystem) = ascii_symbols(typeof(pnsystem))
function ascii_symbols(::Type{T}) where {T<:PNSystem}
    error("`ascii_symbols` is not yet defined for PNSystem subtype `$T`.")
end

"""
    pnsystem::PNSystem(; kwargs...)

State-modifying copy constructor for `PNSystem` objects.

Note that this cannot modify the type's parameters, including the number type `NT`, the
state type `ST`, or the `PNOrder` of the system.  However, it can modify any of the state
variables by symbol or by ASCII symbol.  This function will raise an AssertionError if
any of the keys in `kwargs` is not a valid symbol for the given `PNSystem` type.

```jldoctest
julia> using PostNewtonian: BBH

julia> pnsystem = BBH(ones(14)/2; PNOrder=7//2)
BBH{Vector{Float64}, 7//2}([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

julia> pnsystem2 = pnsystem(M₁=0.2, M₂=0.8, chi1x=0.1)
BBH{Vector{Float64}, 7//2}([0.2, 0.8, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
```
"""
function (pnsystem::PNSystem{NT,PNOrder,ST})(; kwargs...) where {NT,PNOrder,ST}
    all_symbols = Set(symbols(pnsystem)) ∪ Set(ascii_symbols(pnsystem))
    @assert keys(kwargs) ⊆ all_symbols (
        "PNSystem of type $(typeof(pnsystem)) does not have these symbols, which were " *
        "input as keyword arguments:\n    $(setdiff(keys(kwargs), all_symbols))\n" *
        "Maybe you passed `String`s instead of `Symbol`s?\n" *
        "The available symbols for this type are\n" *
        "    $(symbols(pnsystem))\n" *
        "and their ASCII equivalents:\n" *
        "    $(ascii_symbols(pnsystem))"
    )
    state = Tuple(
        get(kwargs, symbol, get(kwargs, ascii_symbol, pnsystem[symbol])) for
        (symbol, ascii_symbol) ∈ zip(symbols(pnsystem), ascii_symbols(pnsystem))
    )
    typeof(pnsystem)(ST(SVector(state)))
end
