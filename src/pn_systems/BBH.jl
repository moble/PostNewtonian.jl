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
struct BBH{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST

    function BBH{NT,ST,PNOrder}(state) where {NT,ST,PNOrder}
        if eachindex(state) != Base.OneTo(14)
            error(
                "The `state` vector for `BBH` must be indexed from 1 to 14; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,ST,PNOrder}(state)
    end
    function BBH(; M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, PNOrder=typemax(Int), kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,ST,PNOrder}(state)
    end
    function BBH(state; Λ₁=0, Λ₂=0, PNOrder=typemax(Int))
        if eachindex(state) != Base.OneTo(14)
            error(
                "The `state` vector for `BBH` must be indexed from 1 to 14; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        @assert Λ₁ == 0
        @assert Λ₂ == 0
        return new{eltype(state),typeof(state),prepare_pn_order(PNOrder)}(state)
    end
end
const BHBH = BBH

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `BBH` systems.
state(pnsystem::BBH) = pnsystem.state
function symbols(::Type{<:BBH})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end
function ascii_symbols(::Type{<:BBH})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end
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

Λ₁(pnsystem::BBH) = zero(pnsystem)
Λ₂(pnsystem::BBH) = zero(pnsystem)

@testitem "PNSystem constructors" begin
    using Quaternionic

    pnA = BBH(;
        M₁=1.0f0, M₂=2.0f0, χ⃗₁=Float32[3.0, 4.0, 5.0], χ⃗₂=Float32[6.0, 7.0, 8.0], v=0.23f0
    )
    @test pnA.state ==
        Float32[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 1.0; 0.0; 0.0; 0.0; 0.23; 0.0]

    pnB = BBH(;
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        v=0.23f0,
        Φ=9.0f0,
    )
    @test pnB.state ==
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
    @test pn1.state ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 0.0]

    pn2 = BBH(;
        M₁=1.0f0,
        M₂=2.0f0,
        χ⃗₁=Float32[3.0, 4.0, 5.0],
        χ⃗₂=Float32[6.0, 7.0, 8.0],
        R=R,
        v=0.23f0,
        Φ=9.0f0,
    )
    @test pn2.state ≈ [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; components(R)...; 0.23; 9.0]

    pn1.state[end] = 9.0f0
    @test pn1.state == pn2.state
end
