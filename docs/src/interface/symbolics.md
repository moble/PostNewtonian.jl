# Symbolic manipulations

It can be useful to evaluate the post-Newtonian expressions with symbolic
arguments.  To do so, we just need to create a `PNSystem` containing a `state`
vector of symbols.  All of the variables defined in
[`PostNewtonian.FundamentalVariables`](@ref "Fundamental variables") and
[`PostNewtonian.DerivedVariables`](@ref "Derived variables") have methods defined
automatically to generate symbols instead of values when called with a symbolic
`PNSystem`.  In turn, any function modified by the
[`@pn_expression`](@ref PostNewtonian.@pn_expression) macro should
also be able to return a symbolic result, including all functions described in
the "[PN expressions](@ref)" section.

For convenience, this package provides [`symbolic_pnsystem`](@ref), which has
all the fundamental variables stored as
[`Symbolics.jl`](https://symbolics.juliasymbolics.org/) variables.  We can
extract the symbols as usual:
```julia
julia> v = PostNewtonian.v(symbolic_pnsystem)
v
```
More generally, we can obtain a symbolic expression for the [binding energy](@ref
PostNewtonian.𝓔) as
```julia
julia> 𝓔(symbolic_pnsystem)
(-1//2)*M*ν*(1 + ((19//8)*ν - (27//8) - (1//24)*(ν^2))*(v^4) + (ν*(χ₁₂ + 6(χₐₗ^2)) + [...]
```
In fact, this is how the [derivative-of-binding-energy function `𝓔′`](@ref
PostNewtonian.𝓔′) is constructed, essentially as
```julia
julia> ∂ᵥ = Differential(v);

julia> expand_derivatives(∂ᵥ(𝓔(symbolic_pnsystem)))
(-1//2)*M*ν*((-90(v^9)*((M₂*λ₁) / M₁ + (M₁*λ₂) / M₂)) / (M^5) + [...]
```

Note that special care is taken to preserve the types of `Irrational`s and the
arguments of certain functions like `sqrt` and `log`.  Ordinarily, Julia will
evaluate these as `Float64`s; to ensure that they remain symbolic, we have to
wrap them in a function that `Symbolics.jl` will know not to bother expanding:
[`PostNewtonian.hold`](@ref).  While manipulating these expressions
symbolically, you'll probably want to leave those `hold` calls as they are.  If
you convert the expressions to code, Julia will compile them away easily, so you
don't need to do anything more.  However, you *can* remove them using
[`PostNewtonian.unhold`](@ref) if the code is in
[`Expr`](https://docs.julialang.org/en/v1/manual/metaprogramming/#Program-representation)
form.
