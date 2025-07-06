# This directory provides the core functionality for the post-Newtonian framework — things
# that we will use to build out the rest of the code.  As such, it is full of little coding
# details that are important for understanding the detailed structure of the code, but not
# for understanding the physics that code is representing.

# These are low-level utilities.
include("utilities/misc.jl")
include("utilities/truncated_series_monoid.jl")
include("utilities/truncated_series_inversion.jl")

# These are types and modules that help to build the post-Newtonian framework.
include("PNTerm.jl")
include("PNExpansion.jl")
include("PNExpressionArithmetic.jl")

# These macros effectively form the interface used directly in the rest of the code.
include("pn_expansion.jl")
include("pn_expression.jl")
include("pn_reference.jl")
