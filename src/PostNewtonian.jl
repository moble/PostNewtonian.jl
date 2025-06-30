module PostNewtonian

using Base: @propagate_inbounds
using FastDifferentiation: Node as FDNode
#using InlineExports: @public, @export  # See below
using Quaternionic: Quaternionic, QuatVec, Rotor, abs2vec, components, normalize, ⋅, ×
using StaticArrays: @MVector, MVector, SVector
using TestItems: @testitem

# While I wait for https://github.com/dalum/InlineExports.jl/pull/2 to be merged, we do the
# following rather than import from the package itself.  Once that is merged, we can add the
# package as a dependency, uncomment the line above, remove this block, and remove the
# `InlineExports.jl` file.
include("core/utilities/InlineExports.jl")
using .InlineExports: @public, @export


include("core/core.jl")
include("pn_systems/pn_systems.jl")
include("literature/literature.jl")
include("pn_expressions/pn_expressions.jl")
include("interface/interface.jl")

end  # module PostNewtonian
