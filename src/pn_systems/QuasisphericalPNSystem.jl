"""
    Quasispherical{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

Abstract type for post-Newtonian systems that are "quasi-spherical" — i.e., where the
eccentricity is negligible, but the system may be precessing (hence on approximately
spherical orbits) while also inspiraling (hence "quasi-").

The built-in subtypes of this type are [`BBH`](@ref), [`BHNS`](@ref), [`NSNS`](@ref), and
[`FDPN`](@ref).  These systems all have an explicit `state` vector containing the
fundamental variables

    M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ

`BHNS`, `NSNS`, and `FDPN` systems also have a [`Λ₂`](@ref) parameter, which is the
quadrupolar tidal coupling parameter of the second object.  `NSNS` and `FDPN` systems also
have a [`Λ₁`](@ref) parameter, which is the quadrupolar tidal coupling parameter of the
first object.

"""
abstract type Quasispherical{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder} end

# The following are methods of functions defined in `state_variables.jl`
state(pnsystem::Quasispherical) = pnsystem.state
function symbols(::Type{<:Quasispherical})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end
function ascii_symbols(::Type{<:Quasispherical})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end
for (i, symbol) ∈ enumerate(symbols(Quasispherical))
    # This will define, e.g., `M₁(pnsystem::Quasispherical) = pnsystem.state[:M₁]`.  We
    # could do this manually, but this is more concise and less error-prone.
    @eval begin
        $(symbol)(pnsystem::Quasispherical) = @inbounds pnsystem.state[i]
        function symbol_index(::Type{T}, ::Val{$symbol}) where {T<:PNSystem}
            index = findfirst(y -> y == $symbol, symbols(T))
            if index === nothing
                error("Type `$(T)` has no symbol `$($symbol)`")
            else
                index
            end
        end
        function Base.getindex(pnsystem::Quasispherical, ::Val{$symbol})
            @inbounds pnsystem.state[i]
        end
        function Base.setindex!(pnsystem::Quasispherical, v, ::Val{$symbol})
            @inbounds pnsystem.state[i] = v
        end
    end
end

χ⃗₁(pnsystem::Quasispherical) = QuatVec(χ⃗₁ˣ(pnsystem), χ⃗₁ʸ(pnsystem), χ⃗₁ᶻ(pnsystem))
χ⃗₂(pnsystem::Quasispherical) = QuatVec(χ⃗₂ˣ(pnsystem), χ⃗₂ʸ(pnsystem), χ⃗₂ᶻ(pnsystem))
R(pnsystem::Quasispherical) = Rotor(Rʷ(pnsystem), Rˣ(pnsystem), Rʸ(pnsystem), Rᶻ(pnsystem))

"""
    BBH{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a binary black hole system.

The `state` vector here holds the fundamental variables `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`,
and `Φ`, with the spins unpacked into three components each, and `R` unpacked into four —
for a total of 14 elements.

The "orbital phase" `Φ` is tracked as the 14th element of the `state` vector.  This is just
the integral of the (scalar) orbital angular frequency `Ω`, and holds little interest for
general systems beyond a convenient description of how "far" the system has evolved.  For
nonprecessing systems, `Φ` would be sufficient to describe the system's position, which is
more completely described by the `Rotor` `R`.
"""
struct BBH{NT,ST,PNOrder} <: Quasispherical{NT,ST,PNOrder}
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
Λ₁(pnsystem::BBH) = zero(pnsystem)
Λ₂(pnsystem::BBH) = zero(pnsystem)

"""
    BHNS{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There is an additional field `Λ₂`
holding the (constant) tidal-coupling parameter of the neutron star.

Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`, and `Λ₂` always
refer to it; `M₁` and `χ⃗₁` always refer to the black hole.  See also [`NSNS`](@ref).
"""
struct BHNS{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST
    Λ₂::NT

    BHNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function BHNS(;
        M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Λ₂, Φ=0, PNOrder=typemax(Int), kwargs...
    )
        NT, ST, PNOrder, state = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₂))
    end
    function BHNS(state; Λ₂, Λ₁=0, PNOrder=typemax(Int))
        if Λ₁ ≠ 0
            error("`BHNS` systems cannot have a nonzero Λ₁; use `NSNS` instead.")
        end
        @assert length(state) == 14
        NT, ST, PNOrder = eltype(state), typeof(state), prepare_pn_order(PNOrder)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₂))
    end
end
Λ₁(pnsystem::BHNS) = zero(pnsystem)
Λ₂(pnsystem::BHNS) = pnsystem.Λ₂

"""
    NSNS{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There are two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
struct NSNS{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST
    Λ₁::NT
    Λ₂::NT

    NSNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function NSNS(;
        M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Λ₁, Λ₂, Φ=0, PNOrder=typemax(Int), kwargs...
    )
        NT, ST, PNOrder, state = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₁), convert(NT, Λ₂))
    end
    function NSNS(state; Λ₁, Λ₂, PNOrder=typemax(Int))
        @assert length(state) == 14
        NT, ST, PNOrder = eltype(state), typeof(state), prepare_pn_order(PNOrder)
        return new{NT,ST,PNOrder}(state, convert(NT, Λ₁), convert(NT, Λ₂))
    end
end
const BNS = NSNS
Λ₁(pnsystem::NSNS) = pnsystem.Λ₁
Λ₂(pnsystem::NSNS) = pnsystem.Λ₂

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
