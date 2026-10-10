"""
    PNSystem{NT, ST, PNOrder}

TODO UPDATE

Base type for all PN systems, such as `BBH`, `BHNS`, and `NSNS`.

These objects encode all essential properties of the binary, including its current state.
As such, they can be used as inputs to the various [fundamental](@ref Fundamental-variables)
and [derived variables](@ref Derived-variables), as well as [PN expressions](@ref) and
[dynamics](@ref Dynamics) functions.

All subtypes should contain a `state` vector holding all of the fundamental variables for
the given type of system.  The parameter `ST` is the type of the `state` vector — for
example, `Vector{Float64}`.  `PNOrder` is a `Rational` giving the order to which PN
expansions should be carried.
"""
abstract type PNSystem{NT,ST<:DenseVector{NT},PNOrder} <: DenseVector{NT} end

abstract type QuasisphericalSystem{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder} end

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

const VecOrPNSystem = Union{AbstractVector,PNSystem}

# Base.eltype(::Type{PNT}) where {NT,PNT<:PNSystem{NT}} = NT
Base.one(::Type{PNT}) where {PNT<:PNSystem} = one(eltype(PNT))
Base.one(x::T) where {T<:PNSystem} = one(T)
Base.zero(::Type{PNT}) where {PNT<:PNSystem} = zero(eltype(PNT))
Base.zero(x::T) where {T<:PNSystem} = zero(T)
Base.float(::Type{PNT}) where {PNT<:PNSystem} = float(eltype(PNT))
Base.float(x::T) where {T<:PNSystem} = float(T)


### Interfaces: https://docs.julialang.org/en/v1/manual/interfaces
# Iteration
Base.iterate(pnsystem::PNSystem) = iterate(state(pnsystem))
Base.iterate(pnsystem::PNSystem, istate) = iterate(state(pnsystem), istate)
Base.IteratorSize(::Type{T}) where {T<:PNSystem} = Base.HasShape{1}()
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
Base.ndims(pnsystem::PNSystem) = ndims(state(pnsystem))
Base.size(pnsystem::PNSystem) = size(state(pnsystem))
Base.size(pnsystem::PNSystem, dim) = size(state(pnsystem), dim)
Base.IteratorEltype(::Type{T}) where {T<:PNSystem} = Base.HasEltype()
Base.eltype(::Type{<:PNSystem{NT}}) where {NT} = NT
Base.isdone(pnsystem::PNSystem) = Base.isdone(state(pnsystem))
Base.isdone(pnsystem::PNSystem, iterstate) = Base.isdone(state(pnsystem), iterstate)

# Indexing

#Base.setindex!(pn::PNSystem, v, i::Int) = Base.@propagate_inbounds setindex!(state(pn), v, i)

Base.@propagate_inbounds Base.getindex(pnsystem::PNSystem, i::Int) = getindex(state(pnsystem), i)
Base.@propagate_inbounds Base.getindex(pnsystem::PNSystem, s::Symbol) = getindex(state(pnsystem), symbol_index(typeof(pnsystem), Val(s)))
Base.@propagate_inbounds Base.setindex!(pnsystem::PNSystem, v, s::Symbol) = setindex!(state(pnsystem), v, symbol_index(typeof(pnsystem), Val(s)))

Base.firstindex(pnsystem::PNSystem) = firstindex(state(pnsystem))
Base.lastindex(pnsystem::PNSystem) = lastindex(state(pnsystem))
Base.eachindex(pnsystem::PNSystem) = eachindex(state(pnsystem))
# Abstract arrays
Base.IndexStyle(::Type{T}) where {T<:PNSystem} = Base.IndexLinear()
#Base.length(pnsystem::PNSystem) = length(state(pnsystem))
# Base.similar(pnsystem::PNSystem) = similar(state(pnsystem))
Base.axes(pnsystem::PNSystem) = axes(state(pnsystem))
# Strided Arrays
Base.strides(pnsystem::PNSystem) = strides(state(pnsystem))
function Base.unsafe_convert(::Type{Ptr{T}}, A::PNSystem) where {T}
    Base.unsafe_convert(Ptr{T}, state(A))
end
Base.elsize(::Type{<:PNSystem{T}}) where {T} = sizeof(T)
Base.stride(pnsystem::PNSystem, k::Int) = stride(state(pnsystem), k)

struct SymbolRange
    start::Symbol
    stop::Symbol
end

Base.:(:)(start::Symbol, stop::Symbol) = SymbolRange(start, stop)

Base.@propagate_inbounds Base.getindex(pnsystem::PNSystem, range::SymbolRange) = getindex(state(pnsystem), symbol_index(typeof(pnsystem), Val(range.start)):symbol_index(typeof(pnsystem), Val(range.stop)))
Base.@propagate_inbounds Base.setindex!(pnsystem::PNSystem, v, range::SymbolRange) = setindex!(state(pnsystem), v, symbol_index(typeof(pnsystem), Val(range.start)):symbol_index(typeof(pnsystem), Val(range.stop)))

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
    causes_domain_error!(u̇, p)

Ensure that these parameters correspond to a physically valid set of PN parameters.

If the parameters are not valid, this function should modify `u̇` to indicate that the
current step is invalid.  This is done by filling `u̇` with `NaN`s, which will be detected
by the ODE solver and cause it to try a different (smaller) step size.

Currently, the only check that is done is to test that these parameters result in a PN
parameter v>0.  In the future, this function may be expanded to include other checks.
"""
function causes_domain_error!(u̇, p::PNSystem{NT}) where {NT}
    if p[:v] ≤ 0  # If this is expanded, document the change in the docstring.
        u̇ .= convert(NT, NaN)
        true
    else
        false
    end
end

"""
    prepare_system
"""
function prepare_system(T::Type{<:PNSystem}; PNOrder=max_pn_order, kwargs...)
    state = pack_state(T; kwargs...)
    ST = typeof(state)
    NT = eltype(ST)
    PNOrder = prepare_pn_order(PNOrder)
    return (NT, ST, PNOrder, state)
end

"""
    prepare_pn_order(PNOrder)

Convert the input to a half-integer of type `Rational{Int}`.

If `PNOrder` is larger than `max_pn_order`, it is set to `max_pn_order`, to avoid overflow
when computing the order index.
"""
function prepare_pn_order(PNOrder)
    if PNOrder < max_pn_order
        round(Int, 2PNOrder)//2
    else
        max_pn_order
    end
end

function StaticArrays.SVector(pnsystem::P) where {P<:PNSystem}
    return SVector{length(pnsystem), eltype(pnsystem)}(pnsystem.state)
end
