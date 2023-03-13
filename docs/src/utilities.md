# Utilities

## Manipulating ODE solutions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/combine_solutions.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```

See also the [Termination criteria](#termination_criteria) page.


## Irrational constants

These quantities are constants that appear in PN expressions, so they are not
exported, but can be used by importing them explicitly or by using the fully
qualified names.  They are defined here as `Irrational`s.  This means that
Julia *can* convert them to float types as necessary.  Unfortunately, by
default Julia converts to `Float64`.  For example, `BigFloat(2ζ3)` will be a
`BigFloat`, but will only have the precision of a `Float64`, because `2ζ3` is
converted first.  To get full precision, you'll need to do things like
`2BigFloat(ζ3)`.

One approach to avoiding this is to explicitly redefine these constants as
floats of the desired precision, using `let` to essentially overwrite the name:
```julia
function foo(x)
    let ζ3=oftype(x, ζ3)
        2ζ3 + x
    end
end
```
Inside the `let` block, `ζ3` is no longer an `Irrational`; it has been converted
to whatever number type `x` is.  Thus, when multiplying by 2, it is not
converted to a `Float64`; its precision matches that of `x`.

This can be quite awkward, so the macro
[`PostNewtonian.@compute_pn_variables`](@ref) is provided to (among other
things) automatically search for all `Irrational`s and replace them with the
appropriate float values.

```@autodocs
Modules = [PostNewtonian]
Filter = t -> typeof(t) !== Irrational{:apery}
Pages   = ["utilities/mathconstants.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
```@docs
PostNewtonian.ζ3
```


## Miscellaneous

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/misc.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
