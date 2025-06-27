"""
    BBH{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a binary black hole system.

The `state` vector here holds the fundamental state variables characterizing the masses,
spins, orientation, velocity, and orbital phase of the system.  The spins unpacked into
three components each.  The orientation is described by the four components of the `Rotor`
`R`.  This gives us a total of 14 elements:

    M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ

The "orbital phase" `Φ` is tracked as the 14th element of the `state` vector.  This is just
the integral of the (scalar) orbital angular frequency `Ω`, and holds little interest for
general systems beyond a convenient description of how "far" the system has evolved.  For
nonprecessing systems, `Φ` would be sufficient to describe the system's position, which is
more completely described by the `Rotor` `R`.  However, for precessing systems, it is
difficult to extract this quantity from `R`.
"""
@export struct BBH{NT,PNOrder,ST} <: Quasispherical{NT,PNOrder,ST}
    state::ST

    function BBH{NT,PNOrder,ST}(state) where {NT,PNOrder,ST}
        if eachindex(state) != Base.OneTo(14)
            error(
                "The `state` vector for `BBH` must be indexed from 1 to 14; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,PNOrder,ST}(state)
    end
    function BBH(; M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, PNOrder=typemax(Int), Λ₁=0, Λ₂=0)
        if Λ₁ != 0 || Λ₂ != 0
            error(
                "`BBH` does not support tidal-coupling parameters `Λ₁` or `Λ₂`; " *
                "use `BHNS` or `NSNS` instead.",
            )
        end
        (NT, PNOrder, state) = prepare_Quasispherical(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,PNOrder,typeof(state)}(state)
    end
    function BBH(state; PNOrder=typemax(Int))
        if eachindex(state) != Base.OneTo(14)
            error(
                "The `state` vector for `BBH` must be indexed from 1 to 14; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        NT, PNOrder, ST = eltype(state), prepare_pn_order(PNOrder), typeof(state)
        return new{NT,PNOrder,ST}(state)
    end
end
@public const BHBH = BBH

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `BBH` systems.
state(pnsystem::BBH) = pnsystem.state
symbols(::Type{<:BBH}) = symbols(Quasispherical)
ascii_symbols(::Type{<:BBH}) = ascii_symbols(Quasispherical)
for (i, symbol) ∈ enumerate(symbols(BBH))
    # This will define, e.g., `M₁(pnsystem::BBH) = pnsystem.state[1]`.  We
    # could do this manually, but this is more concise and less error-prone.
    @eval begin
        $(symbol)(pnsystem::BBH) = @inbounds pnsystem.state[$i]
        function symbol_index(::Type{T}, ::Val{Symbol($symbol)}) where {T<:BBH}
            $i
        end
    end
end

# Define any state-variable methods we may need for variables that are not actually in the
# state vector.
Λ₁(pnsystem::BBH) = zero(pnsystem)
Λ₂(pnsystem::BBH) = zero(pnsystem)

@testitem "BBH constructors" begin
    using PostNewtonian: state
    using Quaternionic

    pnA = BBH(;
        M₁=1.0f0, M₂=2.0f0, χ⃗₁=Float32[3.0, 4.0, 5.0], χ⃗₂=Float32[6.0, 7.0, 8.0], v=0.23f0
    )
    @test eltype(pnA) == Float32
    @test state(pnA) ==
        Float32[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 1.0; 0.0; 0.0; 0.0; 0.23; 0.0]

    pnB = BBH(;
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Φ=9.0f0,
    )
    @test state(pnB) ==
        Float32[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 1.0; 0.0; 0.0; 0.0; 0.23; 9.0]

    R = randn(RotorF32)
    pn1 = BBH(;
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
    )
    @test state(pn1) ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 0.0]

    pn2 = BBH(;
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Φ=9.0f0,
    )
    @test state(pn2) ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 9.0]

    state(pn1)[end] = 9.0f0
    @test state(pn1) == state(pn2)
end
