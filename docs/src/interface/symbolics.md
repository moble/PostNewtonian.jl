# Symbolic manipulations

It can be useful to evaluate the post-Newtonian expressions with
symbolic arguments.  To do so, we just need to create a `PNSystem`
containing a `state` vector of symbols.  All of the variables defined
in [`PostNewtonian.FundamentalVariables`](@ref "Fundamental
variables") and [`PostNewtonian.DerivedVariables`](@ref "Derived
variables") have methods defined automatically to generate symbols
instead of values when called with a symbolic `PNSystem`.  In turn,
any function modified by the [`@pn_expression`](@ref
PostNewtonian.@pn_expression) macro should also be able to return a
symbolic result, including all functions described in the "[PN
expressions](@ref)" section.

For convenience, an extension to this package also provides
[`SymbolicPNSystem`](@ref PostNewtonian.SymbolicPNSystem), which
produces a `PNSystem` with all the fundamental variables stored as
[`Symbolics.jl`](https://symbolics.juliasymbolics.org/) variables — as
long as the `Symbolics` package is also loaded.  We can extract the
symbols as usual:

```julia
julia> using PostNewtonian, Symbolics

julia> symbolic_pnsystem = PostNewtonian.SymbolicPNSystem()
PostNewtonian.SymbolicPNSystem{Vector{Num}, 9223372036854775805//2, Num}(Num[M₁, M₂, χ₁ˣ, χ₁ʸ, χ₁ᶻ, χ₂ˣ, χ₂ʸ, χ₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ], Λ₁, Λ₂)

julia> v = PostNewtonian.v(symbolic_pnsystem)
v
```

More generally, we can obtain a symbolic expression for the [binding
energy](@ref PostNewtonian.𝓔) as

```julia
julia> 𝓔(symbolic_pnsystem)
(-1//2)*(1 + (v^2)*(-(3//4) - (1//12)*ν) + (v^3)*((14//3)*sₗ + 2δ*σₗ) + [...]
```

In fact, this is how the [derivative-of-binding-energy function
`𝓔′`](@ref PostNewtonian.𝓔′) used to be constructed, essentially as

```julia
julia> ∂ᵥ = Differential(v);

julia> expand_derivatives(∂ᵥ(𝓔(symbolic_pnsystem)))
-(1 + (v^2)*(-(3//4) - (1//12)*ν) + (v^3)*((14//3)*sₗ + 2δ*σₗ) + [...]
```

Note that special care is taken to preserve the types of `Irrational`s
and the arguments of certain functions like `sqrt` and `log`.
Ordinarily, Julia will evaluate these as `Float64`s; to ensure that
they remain symbolic, we have to wrap them in a function that
`Symbolics.jl` will know not to bother expanding:
[`PostNewtonian.hold`](@ref).  While manipulating these expressions
symbolically, you'll probably want to leave those `hold` calls as they
are.  If you convert the expressions to code, Julia will compile them
away easily, so you don't need to do anything more.  However, you
*can* remove them using [`PostNewtonian.unhold`](@ref) if the code is
in
[`Expr`](https://docs.julialang.org/en/v1/manual/metaprogramming/#Program-representation)
form.

```@docs
PostNewtonian.hold
PostNewtonian.unhold
PostNewtonian.SymbolicPNSystem
```
