# This directory provides the core functionality for the post-Newtonian framework — things
# that we will use to build out the rest of the code.  As such, it is full of little coding
# details that are important for understanding the detailed structure of the code, but not
# for understanding the physics that code is representing.

include("utilities/misc.jl")
include("utilities/truncated_series_monoid.jl")
include("utilities/truncated_series_inversion.jl")
include("PNTerm.jl")
include("PNExpansion.jl")
include("PNExpressionArithmetic.jl")
include("PNExpression.jl")
include("PNReference.jl")
