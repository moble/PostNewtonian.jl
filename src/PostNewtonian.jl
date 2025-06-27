module PostNewtonian

using Base: @propagate_inbounds
using FastDifferentiation: Node as FDNode
using InlineExports: @public, @export
using Quaternionic: Quaternionic, QuatVec, Rotor, abs2vec, components, normalize, ⋅, ×
using StaticArrays: @MVector, MVector, SVector
using TestItems: @testitem

include("core/core.jl")
include("pn_systems/pn_systems.jl")
include("literature/literature.jl")
include("pn_expressions/pn_expressions.jl")
include("interface/interface.jl")

end  # module PostNewtonian
