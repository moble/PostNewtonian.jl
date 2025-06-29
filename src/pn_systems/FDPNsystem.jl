"""
    FDPNSystem{NT, PN, PNOrder} <: PNSystem{FDNode, Vector{FDNode}, PNOrder}

A `PNSystem` that contains information as variables from
[`FastDifferentiation.jl`](https://docs.juliahub.com/General/FastDifferentiation/stable/).

Note that this type also involves the type parameter `PN`, which is actually the type of a
`PNSystem`, and its type parameter `NT`, which will be the number type of actual numbers
that eventually get fed into (and will be passed out from) functions that use this system.

One important example of what this type is used for is computing the derivative of the
orbital binding energy, `𝓔′` — and in particular, for generating the corresponding function
method to apply to a given `PNSystem`.
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

Base.eltype(::FDPNSystem{FT}) where {FT} = FT

@testitem "FDPNSystem" begin
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
        end
    end
end
