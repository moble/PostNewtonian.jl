function Base.ismutable(pnsystem::PNSystem{NT,PNOrder,ST}) where {NT,PNOrder,ST}
    ismutable(state(pnsystem))
end
function Base.ismutabletype(::Type{<:PNSystem{NT,PNOrder,ST}}) where {NT,PNOrder,ST}
    ismutabletype(ST)
end

### Symbol-based indexing
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
function Base.getindex(pnsystem::T, s::Val{S}) where {T<:PNSystem,S}
    # If `S` is not actually a symbol in `pnsystem`, `symbol_index` will error, so we know
    # that the `index` is inbounds if it returns.
    index = symbol_index(T, s)
    @inbounds state(pnsystem)[index]
end
Base.setindex!(pnsystem::PNSystem, v, s::Symbol) = setindex!(pnsystem, v, Val(s))
function Base.setindex!(pnsystem::T, v, ::Val{S}) where {T<:PNSystem,S}
    index = symbol_index(T, Val(S))
    @inbounds setindex!(state(pnsystem), v, index)
end

### Interfaces: https://docs.julialang.org/en/v1/manual/interfaces
# Iteration
# Base.iterate(pnsystem::PNSystem) = iterate(state(pnsystem))
# Base.iterate(pnsystem::PNSystem, state) = iterate(state(pnsystem), state)
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
@propagate_inbounds Base.getindex(pnsystem::PNSystem, i::T) where {T} = getindex(
    state(pnsystem), i
)
@propagate_inbounds Base.setindex!(pn::PNSystem, v, i::T) where {T} = setindex!(
    state(pn), v, i
)
@propagate_inbounds Base.getindex(pnsystem::PNSystem, i...) = getindex(
    state(pnsystem), i...
)
@propagate_inbounds Base.setindex!(pn::PNSystem, v, i...) = setindex!(state(pn), v, i...)
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
function Base.unsafe_convert(::Type{Ptr{ST}}, pnsystem::PNSystem{ST}) where {ST}
    Base.unsafe_convert(Ptr{ST}, state(pnsystem))
end
Base.elsize(::Type{<:PNSystem{T}}) where {T} = sizeof(T)
Base.stride(pnsystem::PNSystem, k::Int) = stride(state(pnsystem), k)

# NamedTuple interface
function Base.convert(::Type{NamedTuple}, pnsystem::PNSystem{N,P,S}) where {N,P,S}
    NamedTuple{symbols(pnsystem),N}(state(pnsystem))
end
Base.keys(pnsystem::PNSystem) = symbols(pnsystem)
function Base.pairs(pnsystem::PNSystem{N,P,S}) where {N,P,S}
    (s=>v for (s, v) ∈ zip(symbols(pnsystem), state(pnsystem)))
end

Base.@propagate_inbounds function Base.getindex(
    pnsystem::PNSystem, s::AbstractVector{Symbol}
)
    # We trick broadcasting into treating `pnsystem` as a scalar by putting it into a
    # 1-tuple.
    getindex.((pnsystem,), s)
end

# # Allow copying LArray of uninitialized data, as with regular Array
Base.copy(pnsystem::PNSystem) = typeof(pnsystem)(copy(state(pnsystem)))
Base.copyto!(x::PNSystem, y::PNSystem) = copyto!(state(x), state(y))

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

# Broadcasting checks for aliasing with Base.dataids but the fallback
# for AbstractArrays is very slow. Instead, we just call dataids on the
# wrapped state
Base.dataids(pnsystem::PNSystem) = Base.dataids(state(pnsystem))

## Misc

# function ArrayInterface.restructure(
#     x::LArray{T,N,D,Syms}, y::LArray{T2,N2,D2,Syms}
# ) where {T,N,D,T2,N2,D2,Syms}
#     reshape(y, size(x)...)
# end

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

@testitem "Vector interface" begin
    using PostNewtonian: state

    for pnsystem ∈ (BBH(randn(14), 7//2), BHNS(randn(15), 7//2), NSNS(randn(16), 7//2))
        @test_throws ErrorException symbol_index(typeof(pnsystem), Val(:nonexistent_symbol))
        @test symbol_index(typeof(pnsystem), Val(:M₁)) == 1
        @test symbol_index(typeof(pnsystem), Val(:M1)) == 1
        @test pnsystem[:M₁] == pnsystem[:M1] == pnsystem[1] == pnsystem.state[1]
        @test symbol_index(typeof(pnsystem), Val(:M₂)) == 2
        @test symbol_index(typeof(pnsystem), Val(:M2)) == 2
        @test pnsystem[:M₂] == pnsystem[:M2] == pnsystem[2] == pnsystem.state[2]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₁ˣ)) == 3
        @test symbol_index(typeof(pnsystem), Val(:chi1x)) == 3
        @test pnsystem[:χ⃗₁ˣ] == pnsystem[:chi1x] == pnsystem[3] == pnsystem.state[3]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₁ʸ)) == 4
        @test symbol_index(typeof(pnsystem), Val(:chi1y)) == 4
        @test pnsystem[:χ⃗₁ʸ] == pnsystem[:chi1y] == pnsystem[4] == pnsystem.state[4]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₁ᶻ)) == 5
        @test symbol_index(typeof(pnsystem), Val(:chi1z)) == 5
        @test pnsystem[:χ⃗₁ᶻ] == pnsystem[:chi1z] == pnsystem[5] == pnsystem.state[5]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₂ˣ)) == 6
        @test symbol_index(typeof(pnsystem), Val(:chi2x)) == 6
        @test pnsystem[:χ⃗₂ˣ] == pnsystem[:chi2x] == pnsystem[6] == pnsystem.state[6]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₂ʸ)) == 7
        @test symbol_index(typeof(pnsystem), Val(:chi2y)) == 7
        @test pnsystem[:χ⃗₂ʸ] == pnsystem[:chi2y] == pnsystem[7] == pnsystem.state[7]
        @test symbol_index(typeof(pnsystem), Val(:χ⃗₂ᶻ)) == 8
        @test symbol_index(typeof(pnsystem), Val(:chi2z)) == 8
        @test pnsystem[:χ⃗₂ᶻ] == pnsystem[:chi2z] == pnsystem[8] == pnsystem.state[8]
        @test symbol_index(typeof(pnsystem), Val(:Rʷ)) == 9
        @test symbol_index(typeof(pnsystem), Val(:Rw)) == 9
        @test pnsystem[:Rʷ] == pnsystem[:Rw] == pnsystem[9] == pnsystem.state[9]
        @test symbol_index(typeof(pnsystem), Val(:Rˣ)) == 10
        @test symbol_index(typeof(pnsystem), Val(:Rx)) == 10
        @test pnsystem[:Rˣ] == pnsystem[:Rx] == pnsystem[10] == pnsystem.state[10]
        @test symbol_index(typeof(pnsystem), Val(:Rʸ)) == 11
        @test symbol_index(typeof(pnsystem), Val(:Ry)) == 11
        @test pnsystem[:Rʸ] == pnsystem[:Ry] == pnsystem[11] == pnsystem.state[11]
        @test symbol_index(typeof(pnsystem), Val(:Rᶻ)) == 12
        @test symbol_index(typeof(pnsystem), Val(:Rz)) == 12
        @test pnsystem[:Rᶻ] == pnsystem[:Rz] == pnsystem[12] == pnsystem.state[12]
        @test symbol_index(typeof(pnsystem), Val(:v)) == 13
        @test pnsystem[:v] == pnsystem[13] == pnsystem.state[13]
        @test symbol_index(typeof(pnsystem), Val(:Φ)) == 14
        @test symbol_index(typeof(pnsystem), Val(:Phi)) == 14
        @test pnsystem[:Φ] == pnsystem[:Phi] == pnsystem[14] == pnsystem.state[14]
    end
    let
        pnsystem = BHNS(randn(15), 7//2)
        @test symbol_index(typeof(pnsystem), Val(:Λ₂)) == 15
        @test symbol_index(typeof(pnsystem), Val(:Lambda2)) == 15
        # @test pnsystem[:Λ₂] == pnsystem[:Lambda2] == pnsystem[15] == pnsystem.state[15]
    end
    let
        pnsystem = NSNS(randn(16), 7//2)
        @test symbol_index(typeof(pnsystem), Val(:Λ₁)) == 15
        @test symbol_index(typeof(pnsystem), Val(:Lambda1)) == 15
        @test pnsystem[:Λ₁] == pnsystem[:Lambda1] == pnsystem[15] == pnsystem.state[15]
        @test symbol_index(typeof(pnsystem), Val(:Λ₂)) == 16
        @test symbol_index(typeof(pnsystem), Val(:Lambda2)) == 16
        @test pnsystem[:Λ₂] == pnsystem[:Lambda2] == pnsystem[16] == pnsystem.state[16]
    end

    for FT ∈ (Float16, Float32, Float64)
        for PNOrder ∈ (0, 1, 3//2, 2, 5//2, 3, 7//2, 4)
            pnsystem = BBH(randn(FT, 14), PNOrder)
            @test Base.ismutable(pnsystem)
            @test Base.ismutabletype(typeof(pnsystem))
            @test collect(pnsystem) == state(pnsystem)
            @test length(pnsystem) == 14
            @test ndims(pnsystem) == 1
            @test size(pnsystem) == (14,)
            @test size(pnsystem, 1) == 14
            @test axes(pnsystem) == (Base.OneTo(14),)
            @test Base.eltype(pnsystem) == FT
            @test Base.eltype(typeof(pnsystem)) == FT
            @test Base.eltype(state(pnsystem)) == FT
            @test Base.elsize(pnsystem) == sizeof(FT)
            @test Base.elsize(typeof(pnsystem)) == sizeof(FT)
            for (i, v) ∈ enumerate(pnsystem)
                @test v == pnsystem[i] == pnsystem.state[i]
            end
        end
    end
end
