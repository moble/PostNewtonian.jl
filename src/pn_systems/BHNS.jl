"""
    BHNS{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref), with an additional field `Λ₂` holding
the (constant) tidal-coupling parameter of the neutron star.

Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`, and `Λ₂` always
refer to it; `M₁` and `χ⃗₁` always refer to the black hole.  (It's "BHNS", not "NSBH".)  See
also [`NSNS`](@ref).
"""
@export struct BHNS{NT,PNOrder,ST} <: Quasispherical{NT,PNOrder,ST}
    state::ST

    function BHNS{NT,PNOrder,ST}(state) where {NT,PNOrder,ST}
        if eachindex(state) != Base.OneTo(15)
            error(
                "The `state` vector for `BHNS` must be indexed from 1 to 15; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,PNOrder,ST}(state)
    end
    function BHNS(; M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, Λ₁=0, Λ₂, PNOrder=typemax(Int))
        if Λ₁ != 0
            error(
                "`BHNS` does not support a tidal-coupling parameter `Λ₁` for object 1; " *
                "use `NSNS` instead.",
            )
        end
        (NT, PNOrder, state) = prepare_Quasispherical(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        state = vcat(state, Λ₂)
        return new{eltype(state),PNOrder,typeof(state)}(state)
    end
end

# The following are methods of functions defined in `vector_interface.jl` and
# `state_variables.jl`, specialized for `BBH` systems.
Base.length(pnsystem::BHNS) = 15  # Specialize this just for efficiency
state(pnsystem::BHNS) = pnsystem.state
function symbols(::Type{<:BHNS})
    (symbols(Quasispherical)..., :Λ₂)
end
function ascii_symbols(::Type{<:BHNS})
    (ascii_symbols(Quasispherical)..., :Lambda2)
end
for (i, (symbol, ascii_symbol)) ∈ enumerate(zip(symbols(BHNS), ascii_symbols(BHNS)))
    # We could do this manually, but this is more concise and less error-prone.
    @eval begin
        # Define, e.g., `M₁(pnsystem::BHNS) = pnsystem.state[1]`.
        $(symbol)(pnsystem::BHNS) = @inbounds pnsystem.state[$i]

        # Specialize `symbol_index` for Val{:M₁}, Val{:M₂}, etc.
        symbol_index(::Type{T}, ::Val{$(QuoteNode(symbol))}) where {T<:BHNS} = $i
    end
    if symbol ≠ ascii_symbol
        @eval begin
            # Specialize `symbol_index` for Val{:M1}, Val{:M2}, etc.
            symbol_index(::Type{T}, ::Val{$(QuoteNode(ascii_symbol))}) where {T<:BHNS} = $i
        end
    end
end

# Define any state-variable methods we may need for variables that are not actually in the
# state vector.
Λ₁(pnsystem::BHNS) = zero(pnsystem)

@testitem "BHNS constructors" begin
    using PostNewtonian: state
    using Quaternionic

    # minimal constructor: default Φ=0, R=Rotor(1)
    pnA = BHNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Λ₂=4.0f0,
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
        4.0
    ]

    # explicit orbital phase
    pnB = BHNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Φ=9.0f0,
        Λ₂=4.0f0,
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
        4.0
    ]

    # custom rotor, default Φ
    R = randn(RotorF32)
    pn1 = BHNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Λ₂=4.0f0,
    )
    @test state(pn1) ≈
        [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 0.0; 4.0]

    # custom rotor and Φ
    pn2 = BHNS(
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Φ=9.0f0,
        Λ₂=4.0f0,
    )
    @test state(pn2) ≈
        [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 9.0; 4.0]

    # mutating the second-to-last element (Φ) to match pn2
    pn1[end - 1] = 9.0f0
    @test state(pn1) == state(pn2)
end
