# Utilities

## Macros

Some of the most useful features of this package are the macros allowing us to
write PN expressions in fairly natural form, without worrying about calculating
all the variables needed for each expression, or manually accounting for the
various PN orders to which we may need to truncate PN expansions.  To achieve
this, we rely primarily on two macros.

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/macros.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Manipulating ODE solutions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/combine_solutions.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Termination criteria

Hopefully, it should not be necessary to directly use these termination
criteria.  They are used by default in the [`orbital_evolution`](@ref) function.
But certain particularly extreme physical parameters may lead ODEs that are
difficult to integrate — especially if new PN systems or terms are introduced.
These, or similar functions may be helpful examples of
["callbacks"](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
that can be passed to the ODE integrator.

Note that several of these will issue warnings if the evolution has to be
terminated for particularly bad or suspicious reasons, even if the `quiet` flag
is set to `true`.  See the documentation of the `quiet` argument to the
[`orbital_evolution`](@ref) function for an example of how to use `Logging` to
quiet even the warnings.

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/termination_criteria.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Irrational constants

These quantities are constants that appear in PN expressions, so they are not
exported, but can be used by importing them explicitly or by using the fully
qualified names.  They are defined here as `Irrational`s.  This means that Julia
*can* convert them to float types as necessary.  Unfortunately, by default Julia
converts to `Float64`.  For example, `BigFloat(2ζ3)` will be a `BigFloat`, but
will only have the precision of a `Float64`, because `2ζ3` is converted first.
To get full precision, you'll need to do things like `2BigFloat(ζ3)`.

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

This can be quite awkward, so the macro [`PostNewtonian.@pn_expression`](@ref)
is provided to (among other things) automatically search for all `Irrational`s
and replace them with the appropriate float values.

```@docs
PostNewtonian.γₑ
PostNewtonian.ζ3
PostNewtonian.ln2
PostNewtonian.ln3
PostNewtonian.ln5
PostNewtonian.ln³╱₂
PostNewtonian.ln⁵╱₂
```


## Miscellaneous

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/misc.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
