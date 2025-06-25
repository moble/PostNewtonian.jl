"""
    PNSystem{NT, ST, PNOrder}

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
abstract type PNSystem{NT,ST<:DenseVector{NT},PNOrder} <: DenseVector{NT} end

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
    pn_order(pnsystem::PNSystem)

Return the PN order of the given `pnsystem`.

This is a `Rational{Int}` that indicates the order to which the PN expansions should be
carried out when using the given object.
"""
pn_order(::PNSystem{NT,ST,PNOrder}) where {NT,ST,PNOrder} = PNOrder

"""
    order_index(pnsystem::PNSystem)

Return the order index of the given `pnsystem`.

This is defined as the (one-based) index into an iterable of PN terms starting at 0pN, then
0.5pN, etc.  Specifically, this is defined as `1 + Int(2pn_order(pnsystem))`.
"""
order_index(pn::PNSystem) = 1 + Int(2pn_order(pn))

"""
    max_pn_order

The maximum PN order that can be used without overflowing the `Int` type.
"""
const max_pn_order = (typemax(Int) - 2) // 2

"""
    prepare_pn_order(PNOrder)

Convert the input to a half-integer of type `Rational{Int}`.

If `PNOrder` is larger than `max_pn_order`, it is set to `max_pn_order`, to avoid overflow
when computing the order index.
"""
function prepare_pn_order(PNOrder)
    if PNOrder < max_pn_order
        round(Int, 2PNOrder) // 2
    else
        max_pn_order
    end
end

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
    symbol_index(::Type{T}, s::Symbol) where {T<:PNSystem}
    symbol_index(::Type{T}, ::Val{s}) where {T<:PNSystem}

Return the index of the symbol `s` in the state vector of the given `PNSystem` type `T`.

Note that the default implementation is slow; `symbol_index(::Type{T}, ::Val{s})` should be
overridden for every symbol (and ASCII equivalent, if desired) for concrete `PNSystem`
types.
"""
function symbol_index(::Type{T}, ::Val{S}) where {T<:PNSystem,S}
    index = findfirst(y -> y == S, symbols(T))
    if isnothing(index)
        index = findfirst(y -> y == S, ascii_symbols(T))
    end
    if isnothing(index)
        error(
            "Type `$(T)` has no symbol `:$(S)`.\n" *
            "Its symbols are `$(symbols(T))`.\n" *
            "The ASCII equivalents are `$(ascii_symbols(T))`.\n",
        )
    else
        @warn "Please define `PostNewtonian.symbol_index(::Type{$T}, ::Val{$S})`"
        index
    end
end

Base.getindex(pnsystem::PNSystem, s::Symbol) = getindex(pnsystem, Val(s))
function Base.getindex(pnsystem::T, ::Val{S}) where {T<:PNSystem,S}
    # If `S` is not actually a symbol in `pnsystem`, `symbol_index` will error, so we know
    # that the `index` is inbounds if it returns.
    index = symbol_index(T, Val(S))
    @inbounds state(pnsystem)[index]
end

Base.setindex!(pnsystem::PNSystem, v, s::Symbol) = setindex!(pnsystem, v, Val(s))
function Base.setindex!(pnsystem::T, v, ::Val{S}) where {NT,T<:PNSystem{NT},S}
    index = symbol_index(T, Val(S))
    @inbounds setindex!(state(pnsystem), v, index)
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
function (pnsystem::PNSystem)(; kwargs...)
    all_symbols = Set(symbols(pnsystem)) ∪ Set(ascii_symbols(pnsystem))
    @assert keys(kwargs) ⊆ all_symbols (
        "PNSystem of type $(typeof(pnsystem)) does not have these symbols which were input:\n" *
        "    $(setdiff(keys(kwargs), all_symbols))\n" *
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
    typeof(pnsystem)(state)
end
