"""
    BBH{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a quasispherical binary black hole system.

The `state` vector here holds the fundamental variables `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`,
with the spins unpacked into three components each, and `R` unpacked into four — for a total
of 13 elements.

Optionally, `Φ` may also be tracked as the 14th element of the `state` vector.  This is just
the integral of the orbital angular frequency `Ω`, and holds little interest for general
systems beyond a convenient description of how "far" the system has evolved.
"""
struct BBH{NT,ST<:DenseVector{NT},PNOrder} <: QuasisphericalSystem{NT,ST,PNOrder}
    state::ST

    BBH{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)

    function BBH(; PNOrder=max_pn_order, kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(BBH; PNOrder=PNOrder, kwargs...)
        return new{NT,ST,PNOrder}(state)
    end
    function BBH(state; PNOrder=max_pn_order)
        @assert length(state) == length(symbols(BBH))
        return new{eltype(state),typeof(state),prepare_pn_order(PNOrder)}(state)
    end
end
const BHBH = BBH

function pack_state(::Type{<:BBH}; M₁, M₂, χ⃗₁, χ⃗₂, R=Rotor(1), v, Φ=0)
    [M₁; M₂; vec(QuatVec(χ⃗₁)); vec(QuatVec(χ⃗₂)); components(Rotor(R)); v; Φ]
end

state(pnsystem::BBH) = pnsystem.state

function symbols(::Type{<:BBH})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end

function ascii_symbols(::Type{<:BBH})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end

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
