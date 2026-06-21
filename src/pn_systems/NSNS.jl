"""
    NSNS{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref), with two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
struct NSNS{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST

    function NSNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder}
        if eachindex(state) != Base.OneTo(16)
            error(
                "The `state` vector for `NSNS` must be indexed from 1 to 16; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,ST,PNOrder}(state)
    end
    function NSNS(;
        M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, Λ₁, Λ₂, PNOrder=typemax(Int), kwargs...
    )
        NT, ST, PNOrder, state = prepare_system(;
            M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₁, Λ₂, PNOrder
        )
        return new{NT,ST,PNOrder}(state)
    end
    function NSNS(state; PNOrder=typemax(Int))
        if eachindex(state) != Base.OneTo(16)
            error(
                "The `state` vector for `NSNS` must be indexed from 1 to 16; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        NT, ST, PNOrder = eltype(state), typeof(state), prepare_pn_order(PNOrder)
        return new{NT,ST,PNOrder}(state)
    end
end
const BNS = NSNS

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `NSNS` systems.
state(pnsystem::NSNS) = pnsystem.state
function symbols(::Type{<:NSNS})
    (
        :M₁,
        :M₂,
        :χ⃗₁ˣ,
        :χ⃗₁ʸ,
        :χ⃗₁ᶻ,
        :χ⃗₂ˣ,
        :χ⃗₂ʸ,
        :χ⃗₂ᶻ,
        :Rʷ,
        :Rˣ,
        :Rʸ,
        :Rᶻ,
        :v,
        :Φ,
        :Λ₁,
        :Λ₂,
    )
end
function ascii_symbols(::Type{<:NSNS})
    (
        :M1,
        :M2,
        :chi1x,
        :chi1y,
        :chi1z,
        :chi2x,
        :chi2y,
        :chi2z,
        :Rw,
        :Rx,
        :Ry,
        :Rz,
        :v,
        :Phi,
        :Lambda1,
        :Lambda2,
    )
end
for (i, symbol) ∈ enumerate(symbols(NSNS))
    # This will define, e.g., `M₁(pnsystem::NSNS) = pnsystem.state[1]`.  We
    # could do this manually, but this is more concise and less error-prone.
    @eval begin
        $(symbol)(pnsystem::NSNS) = @inbounds pnsystem.state[$i]
        function symbol_index(::Type{T}, ::Val{Symbol($symbol)}) where {T<:NSNS}
            $i
        end
    end
end

Λ₁(pnsystem::NSNS) = @inbounds pnsystem.state[15]
Λ₂(pnsystem::NSNS) = @inbounds pnsystem.state[16]
