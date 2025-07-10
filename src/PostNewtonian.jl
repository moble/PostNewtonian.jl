module PostNewtonian

using Base: @propagate_inbounds
using FastDifferentiation: FastDifferentiation, Node as FDNode
#using InlineExports: @public, @export  # See below
using IrrationalConstants: @irrational
using MacroTools: MacroTools
using Quaternionic:
    Quaternionic, QuatVec, Rotor, 𝐢, 𝐣, 𝐤, abs2vec, absvec, components, normalize, ⋅, ×
using StaticArrays: MVector, SVector
using TestItems: @testitem

# While I wait for https://github.com/dalum/InlineExports.jl/pull/2 to be merged, we do the
# following rather than import from the package itself.  Once that is merged, we can add the
# package as a dependency, uncomment the line above, remove this block, and remove the
# following file.
include("core/utilities/InlineExports.jl")
using .InlineExports: @public, @export

# These are definitions / aliases that are common in PN literature and/or used throughout
# this package.
@public const 𝒾 = im  # Type this as `\scri<tab>`
@public const γₑ = Base.MathConstants.γ  # Distinguish Euler's constant from `γₚₙ = M/r`
public ζ3  # Defined and documented in `core/utilities/misc.jl`

# We will use these types to ensure that precision is preserved in PN expressions.
@public const ExactReal = Union{Integer,Rational,AbstractIrrational}
@public const ExactNumber = Union{ExactReal,Complex{<:ExactReal}}
@public const ExactIntegerBased = Union{Integer,Rational}

include("pn_systems/pn_systems.jl")
include("core/core.jl")
include("literature/literature.jl")
include("pn_expressions/pn_expressions.jl")
include("interface/interface.jl")

end  # module PostNewtonian
