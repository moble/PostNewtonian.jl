"""
    NSNS{NT, PNOrder, ST} <: PNSystem{NT, PNOrder, ST}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref), with two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
@export struct NSNS{NT,PNOrder,ST} <: Quasispherical{NT,PNOrder,ST}
    state::ST

    function NSNS{NT,PNOrder,ST}(state) where {NT,PNOrder,ST}
        if eachindex(state) != Base.OneTo(16)
            error(
                "The `state` vector for `NSNS` must be indexed from 1 to 16; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,PNOrder,ST}(state)
    end
    function NSNS(; M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, Λ₁, Λ₂, PNOrder=typemax(Int))
        (NT, PNOrder, state) = prepare_Quasispherical(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        state = vcat(state, Λ₁, Λ₂)
        return new{eltype(state),PNOrder,typeof(state)}(state)
    end
    # function NSNS(state, PNOrder=typemax(Int))
    #     if eachindex(state) != Base.OneTo(16)
    #         error(
    #             "The `state` vector for `NSNS` must be indexed from 1 to 16; " *
    #             "input is indexed `$(eachindex(state))`.",
    #         )
    #     end
    #     NT, PNOrder, ST = eltype(state), prepare_pn_order(PNOrder), typeof(state)
    #     return new{NT,PNOrder,ST}(state)
    # end
end
@public const BNS = NSNS

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `NSNS` systems.
state(pnsystem::NSNS) = pnsystem.state
function symbols(::Type{<:NSNS})
    (symbols(Quasispherical)..., :Λ₁, :Λ₂)
end
function ascii_symbols(::Type{<:NSNS})
    (ascii_symbols(Quasispherical)..., :Lambda1, :Lambda2)
end
for (i, (symbol, ascii_symbol)) ∈ enumerate(zip(symbols(NSNS), ascii_symbols(NSNS)))
    # We could do this manually, but this is more concise and less error-prone.
    @eval begin
        # Define, e.g., `M₁(pnsystem::NSNS) = pnsystem.state[1]`.
        $(symbol)(pnsystem::NSNS) = @inbounds pnsystem.state[$i]

        # Specialize `symbol_index` for Val{:M₁}, Val{:M₂}, etc.
        symbol_index(::Type{T}, ::Val{$(QuoteNode(symbol))}) where {T<:NSNS} = $i
    end
    if symbol ≠ ascii_symbol
        @eval begin
            # Specialize `symbol_index` for Val{:M1}, Val{:M2}, etc.
            symbol_index(::Type{T}, ::Val{$(QuoteNode(ascii_symbol))}) where {T<:NSNS} = $i
        end
    end
end

@testitem "NSNS constructors" begin
    using PostNewtonian: state
    using Quaternionic

    # minimal constructor: default Φ=0, R=Rotor(1)
    pnA = NSNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Λ₁=5.0f0,
        Λ₂=6.0f0,
    )
    @test eltype(pnA) == Float32
    @test state(pnA) == Float32[
        1.0;
        2.0;
        3.0;
        4.0;
        5.0;
        6.0;
        7.0;
        8.0;
        1.0;
        0.0;
        0.0;
        0.0;
        0.23;
        0.0;
        5.0;
        6.0
    ]

    # explicit orbital phase
    pnB = NSNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Φ=9.0f0,
        Λ₁=5.0f0,
        Λ₂=6.0f0,
    )
    @test state(pnB) == Float32[
        1.0;
        2.0;
        3.0;
        4.0;
        5.0;
        6.0;
        7.0;
        8.0;
        1.0;
        0.0;
        0.0;
        0.0;
        0.23;
        9.0;
        5.0;
        6.0
    ]

    # custom rotor, default Φ
    R = randn(RotorF32)
    pn1 = NSNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Λ₁=5.0f0,
        Λ₂=6.0f0,
    )
    @test state(pn1) ≈
        [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 0.0; 5.0; 6.0]

    # custom rotor and Φ
    pn2 = NSNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Φ=9.0f0,
        Λ₁=5.0f0,
        Λ₂=6.0f0,
    )
    @test state(pn2) ≈
        [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 9.0; 5.0; 6.0]

    # mutating the second-to-last element (Φ) to match pn2
    pn1[end - 2] = 9.0f0
    @test state(pn1) == state(pn2)
end
