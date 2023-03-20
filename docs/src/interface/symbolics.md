# Symbolic manipulations

It can be useful to evaluate the post-Newtonian expressions with symbolic
arguments.  To do so, we just need to create a `PNSystem` containing a `state`
vector of symbols.  All of the variables defined in
[`PostNewtonian.FundamentalVariables`](@ref Fundamental variables) and
[`PostNewtonian.DerivedVariables`](@ref Derived variables) have methods defined
automatically to generate symbols instead of values when called with a symbolic
`PNSystem`.  In turn, any function modified by the
[`@pn_expression`](@ref PostNewtonian.@pn_expression) macro should
also be able to return a symbolic result, including all functions described in
the "[PN expressions](@ref)" section.
