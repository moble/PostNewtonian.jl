# TODO: Write the docstrings for all the general state functions
# TODO: Change to subtyping SVector or MVector
#       - similar_type(::Type{MyPNSystem}, ::Type{NewElType})
#       - similar(::Type{MyPNSystem}, ::Type{NewElType}) for MVector
#       - setindex! for MVector
#       - zero-parameter constructor for MVector

"""
    PNSystem{NT, PNOrder}

Base type for all PN systems, such as `BBH`, `BHNS`, and `NSNS`.

These objects encode all essential properties of the binary, including its current state.
As such, they can be used as inputs to the various [fundamental](@ref Fundamental-variables)
and [derived variables](@ref Derived-variables), as well as [PN expressions](@ref) and
[dynamics](@ref Dynamics) functions.

All subtypes should contain a `state` vector holding all of the fundamental variables for
the given type of system.  The parameter `ST` is the type of the `state` vector — for
example, `Vector{Float64}`.  `PNOrder` is a `Rational` giving the order to which PN
expansions should be carried.
"""
abstract type PNSystem{NT,PNOrder} <: Vector{NT} end

"""
    state(pnsystem::PNSystem)

Return the state vector of `pnsystem`, which is a vector of fundamental variables for the
given PN system.

Note that the built-in `PNSystem` subtypes have a `state` field that is a vector, so this
function will just return that vector.  However, that may not always be true for
user-defined subtypes.
"""
function state(::T) where {T<:PNSystem}
    error("`state` is not yet defined for PNSystem subtype `$T`.")
end
Base.vec(pnsystem::PNSystem) = state(pnsystem)

Base.one(::PNSystem{NT}) where {NT} = one(NT)
Base.zero(::PNSystem{NT}) where {NT} = zero(NT)

"""
    G(pnsystem)

Return Newton's gravitational constant for the given `pnsystem`.

By default, the value is one *with the same number type as `pnsystem`*.  It can be
overridden for subtypes of `PNSystem` that use different units or conventions.

However, note that this function should specialize on the number type of `pnsystem`, rather
than just returning the integer `1`, because there may be expressions with factors such as
`G/3` which will immediately convert to `Float64` if `G` is just `1`, so the result will not
have the expected precision.
"""
G(::PNSystem{NT}) where {NT} = one(NT)

"""
    c(pnsystem)

Return the speed of light for the given `pnsystem`.

By default, the value is one *with the same number type as `pnsystem`*.  It can be
overridden for subtypes of `PNSystem` that use different units or conventions.

However, note that this function should specialize on the number type of `pnsystem`, rather
than just returning the integer `1`, because there may be expressions with factors such as
`c/3` which will immediately convert to `Float64` if `c` is just `1`, so the result will not
have the expected precision.
"""
c(::PNSystem{NT}) where {NT} = one(NT)

"""
    symbols(pnsystem::PNSystem)
    symbols(::Type{<:PNSystem})
    ascii_symbols(pnsystem::PNSystem)
    ascii_symbols(::Type{<:PNSystem})

Return a Tuple of symbols corresponding to the variables tracked by `pnsystem`, in the order
in which they are stored in the `state` vector.
"""
symbols(pnsystem::PNSystem) = symbols(typeof(pnsystem))
function symbols(::Type{T}) where {T<:PNSystem}
    error("`symbols`is not yet defined for PNSystem subtype `$T`.")
end
ascii_symbols(pnsystem::PNSystem) = ascii_symbols(typeof(pnsystem))
function ascii_symbols(::Type{T}) where {T<:PNSystem}
    error("`ascii_symbols` is not yet defined for PNSystem subtype `$T`.")
end

# Define all fundamental variables as functions of a PNSystem, but undefined for the
# abstract type.
for var ∈ (
    :M₁,
    :M₂,
    :χ⃗₁,
    :χ⃗₁ˣ,
    :χ⃗₁ʸ,
    :χ⃗₁ᶻ,
    :χ⃗₂,
    :χ⃗₂ˣ,
    :χ⃗₂ʸ,
    :χ⃗₂ᶻ,
    :R,
    :Rʷ,
    :Rˣ,
    :Rʸ,
    :Rᶻ,
    :v,
    :Φ,
    :Λ₁,
    :Λ₂,
)
    @eval begin
        function $var(::T) where {T<:PNSystem}
            error("$var is not (yet) defined for PNSystem subtype `$T`.")
        end
    end
end

const M1 = M₁
const M2 = M₂
const chi1 = χ⃗₁
const chi1x = χ⃗₁ˣ
const chi1y = χ⃗₁ʸ
const chi1z = χ⃗₁ᶻ
const chi2 = χ⃗₂
const chi2x = χ⃗₂ˣ
const chi2y = χ⃗₂ʸ
const chi2z = χ⃗₂ᶻ
# const R = R
const Rw = Rʷ
const Rx = Rˣ
const Ry = Rʸ
const Rz = Rᶻ
# const v = v
const Phi = Φ
const Lambda1 = Λ₁
const Lambda2 = Λ₂

### Interfaces: https://docs.julialang.org/en/v1/manual/interfaces
# Iteration
Base.iterate(pnsystem::PNSystem) = iterate(state(pnsystem))
Base.iterate(pnsystem::PNSystem, state) = iterate(state(pnsystem), state)
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
Base.IteratorSize(::Type{T}) where {T<:PNSystem} = Base.IteratorSize(Vector)
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
Base.size(pnsystem::PNSystem) = size(state(pnsystem))
Base.size(pnsystem::PNSystem, dim) = size(state(pnsystem), dim)
Base.IteratorEltype(::Type{T}) where {T<:PNSystem} = Base.IteratorEltype(Vector)
Base.eltype(::Type{T}) where {NT,T<:PNSystem{NT}} = NT
Base.isdone(pnsystem::PNSystem) = isdone(state(pnsystem))
Base.isdone(pnsystem::PNSystem, iterstate) = isdone(state(pnsystem), iterstate)
# Indexing
Base.getindex(pnsystem::PNSystem, i) = @propagate_inbounds getindex(state(pnsystem), i)
Base.setindex!(pnsys::PNSystem, v, i) = @propagate_inbounds setindex!(state(pnsys), v, i)
Base.firstindex(pnsystem::PNSystem) = firstindex(state(pnsystem))
Base.lastindex(pnsystem::PNSystem) = lastindex(state(pnsystem))
# Abstract arrays
Base.size(pnsystem::PNSystem) = size(state(pnsystem))
Base.IndexStyle(::Type{T}) where {T<:PNSystem} = Base.IndexStyle(Vector)
Base.length(pnsystem::PNSystem) = length(state(pnsystem))
# Base.similar(pnsystem::PNSystem) = similar(state(pnsystem))
Base.axes(pnsystem::PNSystem) = axes(state(pnsystem))
