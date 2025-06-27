function Base.ismutable(pnsystem::PNSystem{NT,PNOrder,ST}) where {NT,PNOrder,ST}
    ismutable(state(pnsystem))
end
function Base.ismutabletype(::Type{<:PNSystem{NT,PNOrder,ST}}) where {NT,PNOrder,ST}
    ismutabletype(ST)
end

"""
    symbol_index(::Type{T}, s::Symbol) where {T<:PNSystem}
    symbol_index(::Type{T}, ::Val{S}) where {T<:PNSystem,S}

Return the index of the symbol `s` in the state vector of the given `PNSystem` type `T`.

Note that the default implementation is slow; `symbol_index(::Type{T}, ::Val{s})` should be
overridden for every symbol (and ASCII equivalent, if desired) for concrete `PNSystem`
types.
"""
@export function symbol_index(::Type{T}, ::Val{S}) where {T<:PNSystem,S}
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
function Base.setindex!(pnsystem::T, v, ::Val{S}) where {T<:PNSystem,S}
    index = symbol_index(T, Val(S))
    @inbounds setindex!(state(pnsystem), v, index)
end

### Interfaces: https://docs.julialang.org/en/v1/manual/interfaces
# Iteration
Base.iterate(pnsystem::PNSystem) = iterate(state(pnsystem))
Base.iterate(pnsystem::PNSystem, state) = iterate(state(pnsystem), state)
Base.IteratorSize(::Type{T}) where {T<:PNSystem} = Base.HasShape{1}()
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
Base.ndims(pnsystem::PNSystem) = ndims(state(pnsystem))
Base.size(pnsystem::PNSystem) = size(state(pnsystem))
Base.size(pnsystem::PNSystem, dim) = size(state(pnsystem), dim)
Base.IteratorEltype(::Type{T}) where {T<:PNSystem} = Base.HasEltype()
# Base.eltype(::Type{<:PNSystem{NT}}) where {NT} = NT  ## Already defined in `PNSystem.jl`
Base.isdone(pnsystem::PNSystem) = isdone(state(pnsystem))
Base.isdone(pnsystem::PNSystem, iterstate) = isdone(state(pnsystem), iterstate)
# Indexing
@propagate_inbounds Base.getindex(pnsystem::PNSystem, i::Int) = getindex(state(pnsystem), i)
@propagate_inbounds Base.setindex!(pn::PNSystem, v, i::Int) = setindex!(state(pn), v, i)
Base.firstindex(pnsystem::PNSystem) = firstindex(state(pnsystem))
Base.lastindex(pnsystem::PNSystem) = lastindex(state(pnsystem))
Base.eachindex(pnsystem::PNSystem) = eachindex(state(pnsystem))
# Abstract arrays
Base.IndexStyle(::Type{T}) where {T<:PNSystem} = Base.IndexLinear()
# Base.length(pnsystem::PNSystem) = length(state(pnsystem))  ## Already defined above
# Base.similar(pnsystem::PNSystem) = similar(state(pnsystem))
Base.axes(pnsystem::PNSystem) = axes(state(pnsystem))
# Strided Arrays
Base.strides(pnsystem::PNSystem) = strides(state(pnsystem))
function Base.unsafe_convert(::Type{Ptr{T}}, A::PNSystem) where {T}
    Base.unsafe_convert(Ptr{T}, state(A))
end
Base.elsize(::Type{<:PNSystem{T}}) where {T} = sizeof(T)
Base.stride(pnsystem::PNSystem, k::Int) = stride(state(pnsystem), k)

# function PreallocationTools.get_tmp(
#     dc::PreallocationTools.DiffCache, u::LArray{T,N,D,Syms}
# ) where {T<:ForwardDiff.Dual,N,D,Syms}
#     nelem = div(sizeof(T), sizeof(eltype(dc.dual_du))) * length(dc.du)
#     if nelem > length(dc.dual_du)
#         PreallocationTools.enlargedualcache!(dc, nelem)
#     end
#     _x = ArrayInterface.restructure(dc.du, reinterpret(T, view(dc.dual_du, 1:nelem)))
#     LabelledArrays.LArray{T,N,D,Syms}(_x)
# end

# function RecursiveArrayTools.recursive_unitless_eltype(
#     a::Type{LArray{T,N,D,Syms}}
# ) where {T,N,D,Syms}
#     LArray{typeof(one(T)),N,D,Syms}
# end

# #####################################
# # NamedTuple compatibility
# #####################################
# ## SLArray to named tuple
# function Base.convert(::Type{NamedTuple}, x::SLArray{S,T,N,L,Syms}) where {S,T,N,L,Syms}
#     tup = NTuple{length(Syms),T}(x.__x)
#     NamedTuple{Syms,typeof(tup)}(tup)
# end
# Base.keys(x::SLArray{S,T,N,L,Syms}) where {S,T,N,L,Syms} = Syms

# ## pairs iterator
# function Base.pairs(x::LArray{T,N,D,Syms}) where {T,N,D,Syms}
#     # (label => getproperty(x, label) for label in Syms) # not type stable?
#     (Syms[i] => x[i] for i ∈ 1:length(Syms))
# end

# function Base.iterate(x::SLArray, args...)
#     iterate(convert(NamedTuple, x), args...)
# end

# #####################################
# # Array Interface
# #####################################
# function Base.print_array(io::IO, w::WignerMatrix{NT,IT}) where {NT,IT<:Rational}
#     Base.print_array(io, parent(w))
# end

# Base.size(x::LArray) = size(getfield(x, :__x))
# Base.@propagate_inbounds Base.getindex(x::LArray, i...) = getfield(x, :__x)[i...]
# Base.@propagate_inbounds function Base.setindex!(x::LArray, y, i...)
#     getfield(x, :__x)[i...] = y
#     return x
# end

# Base.propertynames(::LArray{T,N,D,Syms}) where {T,N,D,Syms} = Syms
# symnames(::Type{LArray{T,N,D,Syms}}) where {T,N,D,Syms} = Syms

# Base.@propagate_inbounds function Base.getproperty(x::LArray, s::Symbol)
#     if s == :__x
#         return getfield(x, :__x)
#     end
#     return getindex(x, Val(s))
# end

# Base.@propagate_inbounds function Base.setproperty!(x::LArray, s::Symbol, y)
#     if s == :__x
#         return setfield!(x, :__x, y)
#     end
#     setindex!(x, y, Val(s))
# end

# Base.@propagate_inbounds Base.getindex(x::LArray, s::Symbol) = getindex(x, Val(s))
# Base.@propagate_inbounds Base.getindex(x::LArray, s::Val) = __getindex(x, s)
# Base.@propagate_inbounds Base.setindex!(x::LArray, v, s::Symbol) = setindex!(x, v, Val(s))

# @generated function Base.setindex!(x::LArray, y, ::Val{s}) where {s}
#     syms = symnames(x)
#     if syms isa NamedTuple
#         idxs = syms[s]
#         return quote
#             Base.@_propagate_inbounds_meta
#             setindex!(getfield(x, :__x), y, $idxs)
#             return x
#         end
#     else # Tuple
#         idx = findfirst(y -> y == s, symnames(x))
#         return quote
#             Base.@_propagate_inbounds_meta
#             setindex!(getfield(x, :__x), y, $idx)
#             return x
#         end
#     end
# end

# Base.@propagate_inbounds function Base.getindex(x::LArray, s::AbstractArray{Symbol,1})
#     [getindex(x, si) for si ∈ s]
# end

# function Base.similar(
#     x::LArray{T,K,D,Syms}, ::Type{S}, dims::NTuple{N,Int}
# ) where {T,K,D,Syms,S,N}
#     tmp = similar(x.__x, S, dims)
#     LArray{S,N,typeof(tmp),Syms}(tmp)
# end

# function StaticArrays.similar_type(
#     ::Type{SLArray{S,T,N,L,Syms}}, T2, ::Size{S}
# ) where {S,T,N,L,Syms}
#     SLArray{S,T2,N,L,Syms}
# end

# # Allow copying LArray of uninitialized data, as with regular Array
# Base.copy(x::LArray) = typeof(x)(copy(getfield(x, :__x)))
# Base.copyto!(x::LArray, y::LArray) = copyto!(getfield(x, :__x), getfield(y, :__x))

# # enable the usage of LAPACK
# function Base.unsafe_convert(::Type{Ptr{T}}, a::LArray{T,N,D,S}) where {T,N,D,S}
#     Base.unsafe_convert(Ptr{T}, getfield(a, :__x))
# end

# Base.convert(::Type{T}, x) where {T<:LArray} = T(x)
# Base.convert(::Type{T}, x::T) where {T<:LArray} = x
# Base.convert(::Type{<:Array}, x::LArray) = convert(Array, getfield(x, :__x))
# function Base.convert(
#     ::Type{AbstractArray{T,N}}, x::LArray{S,N,<:Any,Syms}
# ) where {T,S,N,Syms}
#     LArray{Syms}(convert(AbstractArray{T,N}, getfield(x, :__x)))
# end
# Base.convert(::Type{AbstractArray{T,N}}, x::LArray{T,N}) where {T,N} = x

# function ArrayInterface.restructure(
#     x::LArray{T,N,D,Syms}, y::LArray{T2,N2,D2,Syms}
# ) where {T,N,D,T2,N2,D2,Syms}
#     reshape(y, size(x)...)
# end

# #####################################
# # Broadcast
# #####################################
# struct LAStyle{T,N,L} <: Broadcast.AbstractArrayStyle{N} end
# LAStyle{T,N,L}(x::Val{i}) where {T,N,L,i} = LAStyle{T,N,L}()
# Base.BroadcastStyle(::Type{LArray{T,N,D,L}}) where {T,N,D,L} = LAStyle{T,N,L}()
# function Base.BroadcastStyle(
#     ::LabelledArrays.LAStyle{T,N,L}, ::LabelledArrays.LAStyle{E,N,L}
# ) where {T,E,N,L}
#     LAStyle{promote_type(T, E),N,L}()
# end

# @generated function labels2axes(::Val{t}) where {t}
#     if t isa NamedTuple && all(x -> x isa Union{Integer,UnitRange}, values(t)) # range labelling
#         (Base.OneTo(maximum(Iterators.flatten(v for v ∈ values(t)))),)
#     elseif t isa NTuple{<:Any,Symbol}
#         axes(t)
#     else
#         error(
#             "$t label isn't supported for broadcasting. Try to formulate it in terms of linear indexing.",
#         )
#     end
# end
# function Base.similar(
#     bc::Broadcast.Broadcasted{LAStyle{T,N,L}}, ::Type{ElType}
# ) where {T,N,L,ElType}
#     tmp = similar(Array{ElType}, axes(bc))
#     if axes(bc) != labels2axes(Val(L))
#         return tmp
#     else
#         return LArray{ElType,N,typeof(tmp),L}(tmp)
#     end
# end

# # Broadcasting checks for aliasing with Base.dataids but the fallback
# # for AbstractArrays is very slow. Instead, we just call dataids on the
# # wrapped buffer
# Base.dataids(pnsystem::PNSystem) = Base.dataids(state(pnsystem))

# Base.elsize(::Type{<:LArray{T}}) where {T} = sizeof(T)
