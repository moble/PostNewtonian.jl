"""
    PNSystem{T, PNOrder}

Base type for all PN systems, such as `BBH`, `BHNS`, and `NSNS`.

These objects encode all essential properties of the binary, including its current state.
As such, they can be used as inputs to the various [fundamental](@ref Fundamental-variables)
and [derived variables](@ref Derived-variables), as well as [PN expressions](@ref) and
[dynamics](@ref Dynamics) functions.

All subtypes should contain a `state` vector holding all of the fundamental variables for
the given type of system.  The parameter `T` is the type of the `state` vector — for
example, `Vector{Float64}`.  `PNOrder` is a `Rational` giving the order to which PN
expansions should be carried.
"""
abstract type PNSystem{ST, PNOrder} end

const VecOrPNSystem = Union{AbstractVector, PNSystem}

const pnsystem_symbols = (
    :M₁, :M₂, :χ⃗₁ˣ, :χ⃗₁ʸ, :χ⃗₁ᶻ, :χ⃗₂ˣ, :χ⃗₂ʸ, :χ⃗₂ᶻ, :Rʷ, :Rˣ, :Rʸ, :Rᶻ, :v, :Φ
)

for (i, s) ∈ enumerate(pnsystem_symbols)
    sindex = Symbol("$(s)index")
    @eval const $sindex = $i
end

const χ⃗₁indices = χ⃗₁ˣindex:χ⃗₁ᶻindex
const χ⃗₂indices = χ⃗₂ˣindex:χ⃗₂ᶻindex
const Rindices = Rʷindex:Rᶻindex


Base.eltype(::PNSystem{ST}) where {ST} = eltype(ST)
pn_order(::PNSystem{ST, PNOrder}) where {ST, PNOrder} = PNOrder

order_index(pn::PNSystem) = 1 + Int(2pn_order(pn))

causes_domain_error!(u̇, ::PNSystem{VT}) where {VT<:Vector{Symbolics.Num}} = false
function causes_domain_error!(u̇, p::PNSystem{VT}) where {VT}
    if p.state[vindex] ≤ 0
        u̇ .= convert(eltype(VT), NaN)
        true
    else
        false
    end
end

function prepare_system(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=0,
    PNOrder=typemax(Int)
)
    state = [M₁; M₂; vec(QuatVec(χ⃗₁)); vec(QuatVec(χ⃗₂)); components(Rotor(R)); v; Φ]
    ST = typeof(state)
    PNOrder = prepare_pn_order(PNOrder)
    (ST, PNOrder, state)
end

function prepare_pn_order(PNOrder)
    if PNOrder!=typemax(Int)
        round(Int, 2PNOrder) // 2
    else
        (typemax(Int)-2) // 2
    end
end


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
struct BBH{T, PNOrder} <: PNSystem{T, PNOrder}
    state::T
end
function BBH(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ=0,
    PNOrder=typemax(Int), kwargs...
)
    (T, PNOrder, state) = prepare_system(;M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
    BBH{T, PNOrder}(state)
end
function BBH(state; Λ₁=0, Λ₂=0, PNOrder=typemax(Int))
    @assert length(state) == 14
    BBH{typeof(state), prepare_pn_order(PNOrder)}(state)
end


"""
    BHNS{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a black-hole—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There is an additional field `Λ₂`
holding the (constant) tidal-coupling parameter of the neutron star.

Note that the neutron star is *always* object 2 — meaning that `M₂`, `χ⃗₂`, and `Λ₂` always
refer to it; `M₁` and `χ⃗₁` always refer to the black hole.  See also [`NSNS`](@ref).
"""
struct BHNS{T, PNOrder} <: PNSystem{T, PNOrder}
    state::T
    Λ₂::eltype(T)
end
function BHNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Λ₂, Φ=0,
    PNOrder=typemax(Int), kwargs...
)
    (T, PNOrder, state) = prepare_system(;M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
    BHNS{T, PNOrder}(state, Λ₂)
end
function BHNS(state; Λ₂, Λ₁=0, PNOrder=typemax(Int))
    @assert length(state) == 14
    BHNS{typeof(state), prepare_pn_order(PNOrder)}(state, Λ₂)
end


"""
    NSNS{T, PNOrder}

The [`PNSystem`](@ref) subtype describing a neutron-star—neutron-star binary system.

The `state` vector is the same as for a [`BBH`](@ref).  There are two additional fields `Λ₁`
and `Λ₂` holding the (constant) tidal-coupling parameters of the neutron stars.  See also
[`BHNS`](@ref).
"""
struct NSNS{T, PNOrder} <: PNSystem{T, PNOrder}
    state::T
    Λ₁::eltype(T)
    Λ₂::eltype(T)
end
function NSNS(;
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Λ₁, Λ₂, Φ=0,
    PNOrder=typemax(Int), kwargs...
)
    (T, PNOrder, state) = prepare_system(;M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, PNOrder)
    NSNS{T, PNOrder}(state, Λ₁, Λ₂)
end
function NSNS(state; Λ₁, Λ₂, PNOrder=typemax(Int))
    @assert length(state) == 14
    NSNS{typeof(state), prepare_pn_order(PNOrder)}(state, Λ₁, Λ₂)
end


"""
    SymbolicPNSystem{T, PNOrder}(state, Λ₁, Λ₂)

A `PNSystem` that contains information as variables from
[`Symbolics.jl`](https://symbolics.juliasymbolics.org/).

See also [`symbolic_pnsystem`](@ref) for a particular general instance of this type.
"""
struct SymbolicPNSystem{T, PNOrder} <: PNSystem{T, PNOrder}
    state::T
    Λ₁::eltype(T)
    Λ₂::eltype(T)
end
function SymbolicPNSystem(PNOrder=typemax(Int))
    @variables M₁ M₂ χ⃗₁ˣ χ⃗₁ʸ χ⃗₁ᶻ χ⃗₂ˣ χ⃗₂ʸ χ⃗₂ᶻ Rʷ Rˣ Rʸ Rᶻ v Φ Λ₁ Λ₂
    SymbolicPNSystem{Vector{typeof(M₁)}, prepare_pn_order(PNOrder)}(
        [M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ],
        Λ₁, Λ₂
    )
end

"""
    symbolic_pnsystem

A symbolic `PNSystem` that contains symbolic information for all types of `PNSystem`s.

In particular, note that this object has an (essentially) infinite `PNOrder`, uses the
`TaylorT1` approximant, and has nonzero values for quantities like `Λ₁` and `Λ₂`.  If you
want different choices, you may need to call [`SymbolicPNSystem`](@ref) yourself, or even
construct a different specialized subtype of `PNSystem` (it's not hard).

# Examples
```jldoctest
julia> using PostNewtonian: M₁, M₂, χ⃗₁, χ⃗₂

julia> M₁(symbolic_pnsystem), M₂(symbolic_pnsystem)
(M₁, M₂)

julia> χ⃗₁(symbolic_pnsystem)
χ⃗₁

julia> χ⃗₂(symbolic_pnsystem)
χ⃗₂
```
"""
const symbolic_pnsystem = SymbolicPNSystem()
