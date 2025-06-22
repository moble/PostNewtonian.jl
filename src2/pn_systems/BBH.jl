abstract type AbstractBBH{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder} end

"""
    BBH{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a binary black hole system.

The `state` vector here holds the fundamental variables `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`,
with the spins unpacked into three components each, and `R` unpacked into four — for a total
of 13 elements.

Optionally, `Φ` may also be tracked as the 14th element of the `state` vector.  This is just
the integral of the orbital angular frequency `Ω`, and holds little interest for general
systems beyond a convenient description of how "far" the system has evolved.
"""
struct BBH{NT,ST,PNOrder} <: AbstractBBH{NT,ST,PNOrder}
    state::ST

    BBH{NT,ST,PNOrder}(state) where {NT,ST,PNOrder} = new{NT,ST,PNOrder}(state)
    function BBH(; M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, PNOrder=typemax(Int), kwargs...)
        (NT, ST, PNOrder, state) = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
        return new{NT,ST,PNOrder}(state)
    end
    function BBH(state; Λ₁=0, Λ₂=0, PNOrder=typemax(Int))
        @assert length(state) == 14
        @assert Λ₁ == 0
        @assert Λ₂ == 0
        return new{eltype(state),typeof(state),prepare_pn_order(PNOrder)}(state)
    end
end
const BHBH = BBH

state(pnsystem::BBH) = pnsystem.state

function symbols(::Type{BBH})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ)
end
function ascii_symbols(::Type{BBH})
    (:M1, :M2, :chi1x, :chi1y, :chi1z, :chi2x, :chi2y, :chi2z, :Rw, :Rx, :Ry, :Rz, :v, :Phi)
end

M₁(pnsystem::BBH) = @inbounds state(pnsystem)[1]
M₂(pnsystem::BBH) = @inbounds state(pnsystem)[2]
χ⃗₁(pnsystem::BBH) = @inbounds state(pnsystem)[3:5]
χ⃗₁ˣ(pnsystem::BBH) = @inbounds state(pnsystem)[3]
χ⃗₁ʸ(pnsystem::BBH) = @inbounds state(pnsystem)[4]
χ⃗₁ᶻ(pnsystem::BBH) = @inbounds state(pnsystem)[5]
χ⃗₂(pnsystem::BBH) = @inbounds state(pnsystem)[6:8]
χ⃗₂ˣ(pnsystem::BBH) = @inbounds state(pnsystem)[6]
χ⃗₂ʸ(pnsystem::BBH) = @inbounds state(pnsystem)[7]
χ⃗₂ᶻ(pnsystem::BBH) = @inbounds state(pnsystem)[8]
R(pnsystem::BBH) = @inbounds state(pnsystem)[9:12]
Rʷ(pnsystem::BBH) = @inbounds state(pnsystem)[9]
Rˣ(pnsystem::BBH) = @inbounds state(pnsystem)[10]
Rʸ(pnsystem::BBH) = @inbounds state(pnsystem)[11]
Rᶻ(pnsystem::BBH) = @inbounds state(pnsystem)[12]
v(pnsystem::BBH) = @inbounds state(pnsystem)[13]
Φ(pnsystem::BBH) = @inbounds state(pnsystem)[14]

# Λ₁(pnsystem::BBH) = @inbounds state(pnsystem)[15]
# Λ₂(pnsystem::BBH) = @inbounds state(pnsystem)[16]
