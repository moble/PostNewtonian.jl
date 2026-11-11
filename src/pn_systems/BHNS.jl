"""
    BHNS{NT, ST, PNOrder} <: PNSystem{NT, ST, PNOrder}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref), with an additional field `Λ₂` holding
the (constant) tidal-coupling parameter of the neutron star.

Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`, and `Λ₂` always
refer to it; `M₁` and `χ⃗₁` always refer to the black hole.  (It's "BHNS", not "NSBH".)  See
also [`NSNS`](@ref).
"""
struct BHNS{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST

    function BHNS{NT,ST,PNOrder}(state) where {NT,ST,PNOrder}
        if eachindex(state) != Base.OneTo(15)
            error(
                "The `state` vector for `BHNS` must be indexed from 1 to 15; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        new{NT,ST,PNOrder}(state)
    end
    function BHNS(;
        M₁, M₂, χ⃗₁, χ⃗₂, v, R=Rotor(1), Φ=0, Λ₂, PNOrder=typemax(Int), kwargs...
    )
        NT, ST, PNOrder, state = prepare_system(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₂, PNOrder)
        return new{NT,ST,PNOrder}(state)
    end
    function BHNS(state; PNOrder=typemax(Int))
        if eachindex(state) != Base.OneTo(15)
            error(
                "The `state` vector for `BHNS` must be indexed from 1 to 15; " *
                "input is indexed `$(eachindex(state))`.",
            )
        end
        NT, ST, PNOrder = eltype(state), typeof(state), prepare_pn_order(PNOrder)
        return new{NT,ST,PNOrder}(state)
    end
end

# The following are methods of functions defined in `state_variables.jl`, specialized for
# `BHNS` systems.
state(pnsystem::BHNS) = pnsystem.state
function symbols(::Type{<:BHNS})
    (:M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ, :Λ₂)
end
function ascii_symbols(::Type{<:BHNS})
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
        :Lambda2,
    )
end
for (i, symbol) ∈ enumerate(symbols(BHNS))
    # This will define, e.g., `M₁(pnsystem::BHNS) = pnsystem.state[1]`.  We
    # could do this manually, but this is more concise and less error-prone.
    @eval begin
        $(symbol)(pnsystem::BHNS) = @inbounds pnsystem.state[$i]
        function symbol_index(::Type{T}, ::Val{Symbol($symbol)}) where {T<:BHNS}
            $i
        end
    end
end

Λ₁(pnsystem::BHNS) = zero(pnsystem)
Λ₂(pnsystem::BHNS) = @inbounds pnsystem.state[15]
