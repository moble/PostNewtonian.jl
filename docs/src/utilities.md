# Mass-parameter conversions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["masses.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


# Spin-parameter conversions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["spins.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


# Orbital elements

```@docs
n̂
λ̂
ℓ̂
v
Ω
```


# ODE solutions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["combine_solutions.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


# Irrational constants

These quantities are constants that appear in PN expressions, so they are not
exported, but can be used by importing them explicitly or by using the fully
qualified names.  They are defined here as `Irrational`s.  This means that
Julia *can* convert them to float types as necessary.  Unfortunately, by
default Julia converts to `Float64`.  For example, `BigFloat(2ζ3)` will be a
`BigFloat`, but will only have the precision of a `Float64`, because `2ζ3` is
converted first.  To get full precision, you'll need to do things like
`2BigFloat(ζ3)`.

```@autodocs
Modules = [PostNewtonian]
Filter = t -> typeof(t) !== Irrational{:apery}
Pages   = ["constants.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
```@docs
PostNewtonian.ζ3
```
