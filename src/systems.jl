"""
    PNSystem{T, PNOrder, Expansion}

Base type for all PN systems, such as `BBH`, `BHNS`, and `NSNS`.

The parameter `T` is the basic float type of all variables — usually just `Float64`.
`PNOrder` is a `Rational` giving the order to which PN expansions should be carried.  And
`Expansion` describes the type of expansion — one of `TaylorT1`, `TaylorT4`, or `TaylorT5`.

All subtypes should contain a `state` vector holding all of the fundamental variables for
the given type of system.
"""
abstract type PNSystem{T, PNOrder, Expansion} end

Base.eltype(::PNSystem{T}) where {T} = T
pn_order(::PNSystem{T, PNOrder}) where {T, PNOrder} = PNOrder
expansion_type(::PNSystem{T, P, Expansion}) where {T, P, Expansion} = Expansion

order_index(pn::PNSystem) = 1 + Int(2pn_order(pn))

function prepare_system(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1
)
    state = [M₁; M₂; vec(χ⃗₁); vec(χ⃗₂); components(R); v]
    if !isnothing(Φ)
        state = [state; Φ]
    end
    T = eltype(state)
    PNOrder = if PNOrder!=typemax(Int)
        round(Int, 2PNOrder) // 2
    else
        typemax(Int) // 2
    end
    (T, PNOrder, Expansion, state)
end

"""
    BBH{T, PNOrder, Expansion}

The [`PNSystem`](@ref) subtype describing a binary black hole system.

The `state` vector here holds the fundamental variables `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`,
with the spins unpacked into three components each, and `R` unpacked into four — for a total
of 13 elements.

Optionally, `Φ` may also be tracked as the 14th element of the `state` vector.  This is just
the integral of the orbital angular frequency `Ω`, and holds little interest for general
systems beyond a convenient description of how "far" the system has evolved.
"""
struct BBH{T, PNOrder, Expansion} <: PNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
end
function BBH(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1, kwargs...
)
    (T, PNOrder, Expansion, state) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    BBH{T, PNOrder, Expansion}(state)
end


"""
    BHNS{T, PNOrder, Expansion}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There is an additional field `λ₂`
holding the (constant) Love number, which describes the tidal-coupling parameter of the
neutron star.  Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`,
and `λ₂` always refer to it; `M₁` and `χ⃗₁` always refer to the black hole.
"""
struct BHNS{T, PNOrder, Expansion} <: PNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
    λ₂::T
end
function BHNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, λ₂, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1, kwargs...
)
    (T, PNOrder, Expansion, state) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    BHNS{T, PNOrder, Expansion}(state, λ₂)
end


"""
    NSNS{T, PNOrder, Expansion}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There are two additional fields `λ₁`
and `λ₂` holding the (constant) Love numbers, which describes the tidal-coupling parameter
of the neutron stars.
"""
struct NSNS{T, PNOrder, Expansion} <: PNSystem{T, PNOrder, Expansion}
    state::AbstractVector{T}
    λ₁::T
    λ₂::T
end
function NSNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, λ₁, λ₂, Φ=nothing,
    PNOrder=typemax(Int), Expansion=:TaylorT1, kwargs...
)
    (T, PNOrder, Expansion, state) = prepare_system(;
        M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=nothing,
        PNOrder, Expansion
    )
    NSNS{T, PNOrder, Expansion}(state, λ₁, λ₂)
end

"""
    symbolic_pnsystem(sys)

Construct a symbolic `PNSystem`, specifically of the subtype given by `sys`.


# Examples
```jldoctest
julia> using PostNewtonian: M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, λ₁, λ₂

julia> pn = symbolic_pnsystem(BBH);

julia> M₁(pn), M₂(pn)
(M₁, M₂)

julia> χ⃗₁(pn), χ⃗₂(pn)
(χ⃗₁, χ⃗₂)
```
"""
function symbolic_pnsystem(sys; PNOrder=typemax(Int), Expansion=:TaylorT1)
    @variables M₁ M₂ χ⃗₁ˣ χ⃗₁ʸ χ⃗₁ᶻ χ⃗₂ˣ χ⃗₂ʸ χ⃗₂ᶻ Rʷ Rˣ Rʸ Rᶻ v Φ λ₁ λ₂
    χ⃗₁ = QuatVec(χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ)
    χ⃗₂ = QuatVec(χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ)
    R = Rotor(Rʷ, Rˣ, Rʸ, Rᶻ)
    sys(; M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, λ₁, λ₂)
end
