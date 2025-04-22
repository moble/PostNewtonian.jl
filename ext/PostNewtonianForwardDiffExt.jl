module PostNewtonianForwardDiffExt
isdefined(Base, :get_extension) ? (using ForwardDiff: ForwardDiff) : (import ..ForwardDiff)

import ForwardDiff: Dual
import PostNewtonian: type_converter, FDPNSystem, PNSystem, 𝓔′, 𝓔′code
import FastDifferentiation: FastDifferentiation, Node
import StaticArrays: SVector, MVector
import MacroTools


# These three definitions allow us to call 𝓔′ with ForwardDiff.Dual numbers
function type_converter(::FDPNSystem{Dual{T, V, N}}, x) where {T, V, N}
    return x
end
function type_converter(::FDPNSystem{Dual{T, V, N}}, x::Integer) where {T, V, N}
    return convert(V, x)
end
function type_converter(::FDPNSystem{Dual{T, V, N}}, x::Rational) where {T, V, N}
    return convert(V, x)
end
function type_converter(::FDPNSystem{Dual{T, V, N}}, x::AbstractIrrational) where {T, V, N}
    return convert(V, x)
end

@generated function 𝓔′(
    pnsystem::PNSystem{<:AbstractVector{Dual{T,V,N}}};
    pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum)
) where {T, V, N, PNExpansionReducer}
    𝓔′code(pnsystem, pn_expansion_reducer, Dual{T,V,N}, V)
end


end #module
