"""
    FDPNSystem{NT, PN, PNOrder} <: PNSystem{Node, Vector{Node}, PNOrder}

A `PNSystem` that contains information as variables of type `Node` from
[`FastDifferentiation.jl`](https://docs.juliahub.com/General/FastDifferentiation/stable/).

Note that this type also involves the type parameter `PN`, which is actually the type of a
`PNSystem`, and its type parameter `NT`, which will be the number type of actual numbers
that eventually get fed into (and will be passed out from) functions that use this system.

One important example of what this type is used for is computing the derivative of the
orbital binding energy, `𝓔′` — and in particular, for generating the corresponding function
method to apply to a given `PNSystem`.

!!! warning

    Because of the structure of the type parameters, most methods defined for a general
    `PNSystem` will use `Node` as the number type.  This may be correct in many cases, but
    when `Node` is applied to an integer, for example, it will generally result in a
    `Float64` value, rather than the `NT` type parameter.  For this reason many methods
    will need to be specialized for `FDPNSystem` types, using almost the same definition,
    but with the additional `FD` prefix inserted throughout.
"""
@export struct FDPNSystem{NT,PNOrder,PN<:PNSystem{NT,PNOrder}} <:
               PNSystem{FDNode,PNOrder,Vector{FDNode}}
    state::Vector{FDNode}

    FDPNSystem(pnsystem::PNSystem) = FDPNSystem(typeof(pnsystem))
    function FDPNSystem(::Type{PN}) where {NT,PNOrder,PN<:PNSystem{NT,PNOrder}}
        return new{NT,PNOrder,PN}([FDNode(s) for s ∈ symbols(PN)])
    end
end

state(pnsystem::FDPNSystem) = pnsystem.state

# This is an example of where we need to specialize the method for `FDPNSystem`, as
# explained in the docstring.
Base.eltype(::FDPNSystem{FT}) where {FT} = FT

symbols(pnsystem::FDPNSystem{NT,PNOrder,PN}) where {NT,PNOrder,PN} = symbols(PN)
symbols(::Type{T}) where {NT,PNOrder,PN,T<:FDPNSystem{NT,PNOrder,PN}} = symbols(PN)

function symbol_index(pnsystem::FDPNSystem{NT,PNOrder,PN}, s::Symbol) where {NT,PNOrder,PN}
    symbol_index(PN, Val(s))
end
function symbol_index(::Type{T}, ::Val{S}) where {T<:FDPNSystem,S}
    index = findfirst(y -> y == S, symbols(T))
    if isnothing(index)
        index = findfirst(y -> y == S, ascii_symbols(T))
    end
    if isnothing(index)
        error(
            "Type `$(T)` has no symbol `:$(S)`.\n" *
            "This type's symbols are `$(symbols(T))`.\n" *
            "The ASCII equivalents are `$(ascii_symbols(T))`.\n",
        )
    else
        index
    end
end

# These are more examples of where we need to specialize the method for `FDPNSystem`, as
# explained in the docstring.
constant_convert(::T, x::ExactNumber) where {NT,T<:FDPNSystem{NT}} = NT(x)
constant_convert(::T, x::NT) where {NT,T<:FDPNSystem{NT}} = x

@testitem "FDPNSystem" begin
    using PostNewtonian: constant_convert
    @testset "BBH" begin
        PNOrder = 7//2
        for NT ∈ (Float16, Float64)
            bbh = BBH(randn(NT, 14), PNOrder)
            fdpnsystem = FDPNSystem(bbh)
            @test fdpnsystem isa FDPNSystem{eltype(bbh),PNOrder,typeof(bbh)}
            @test pn_order(fdpnsystem) == PNOrder
            @test eltype(fdpnsystem) == NT
            @test symbols(fdpnsystem) == symbols(bbh)
            @test length(fdpnsystem) == 14
            @test constant_convert(fdpnsystem, π) isa NT
            @test constant_convert(fdpnsystem, π) == NT(π)
        end
    end
    @testset "BHNS" begin
        PNOrder = typemax(Int)
        for NT ∈ (Float16, Float64)
            bhns = BHNS(randn(NT, 15), PNOrder)
            fdpnsystem = FDPNSystem(bhns)
            @test fdpnsystem isa FDPNSystem{eltype(bhns),max_pn_order,typeof(bhns)}
            @test pn_order(fdpnsystem) == max_pn_order
            @test eltype(fdpnsystem) == NT
            @test symbols(fdpnsystem) == symbols(bhns)
            @test length(fdpnsystem) == 15
            @test constant_convert(fdpnsystem, π) isa NT
            @test constant_convert(fdpnsystem, π) == NT(π)
        end
    end
    @testset "NSNS" begin
        PNOrder = 3.5
        for NT ∈ (Float16, Float64)
            nsns = NSNS(randn(NT, 16), PNOrder)
            fdpnsystem = FDPNSystem(nsns)
            @test fdpnsystem isa FDPNSystem{eltype(nsns),7//2,typeof(nsns)}
            @test pn_order(fdpnsystem) == 7//2
            @test eltype(fdpnsystem) == NT
            @test symbols(fdpnsystem) == symbols(nsns)
            @test length(fdpnsystem) == 16
            @test constant_convert(fdpnsystem, π) isa NT
            @test constant_convert(fdpnsystem, π) == NT(π)
        end
    end
end
