module PostNewtonianForwardDiffExt
isdefined(Base, :get_extension) ? (using ForwardDiff: ForwardDiff) : (import ..ForwardDiff)

import ForwardDiff: Dual, valtype
import PostNewtonian: type_converter, PNSystem, PNExpansion, 𝓔′, 𝓔′code
import FastDifferentiation: FastDifferentiation, Node
import StaticArrays: StaticArrays, SVector, MVector
using MacroTools: MacroTools

Base.one(::Type{PNT}) where {PNT<:PNSystem{<:Dual}} = one(valtype(eltype(PNT)))
Base.zero(::Type{PNT}) where {PNT<:PNSystem{<:Dual}} = zero(valtype(eltype(PNT)))
Base.float(::Type{PNT}) where {PNT<:PNSystem{<:Dual}} = float(valtype(eltype(PNT)))

# These three definitions allow us to call 𝓔′ with ForwardDiff.Dual numbers
function type_converter(::PNSystem{Dual{T,V,N}}, x) where {T,V,N}
    return x
end
function type_converter(::PNSystem{Dual{T,V,N}}, x::Integer) where {T,V,N}
    return convert(V, x)
end
function type_converter(::PNSystem{Dual{T,V,N}}, x::Rational) where {T,V,N}
    return convert(V, x)
end
function type_converter(::PNSystem{Dual{T,V,N}}, x::AbstractIrrational) where {T,V,N}
    return convert(V, x)
end

@generated function 𝓔′(
    pnsystem::PNSystem{<:AbstractVector{Dual{T,V,N}}};
    pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum),
) where {T,V,N,PNExpansionReducer}
    𝓔′code(pnsystem, pn_expansion_reducer, Dual{T,V,N}, V)
end

end #module
