# Utilities

## PN terms and PN expansions

The basic building blocks of post-Newtonian theory are the terms and
expansions.  These are used to build up the various expressions that
describe the dynamics of the system.  The terms are the individual
parts of the expansions, while the expansions are the full expressions
that are expanded in powers of ``1/c``.

```@docs
PostNewtonian.PNTerm
PostNewtonian.PNExpansionParameter
PostNewtonian.PNExpansion
```

## Macros

Some of the most useful features of this package are the macros allowing us to
write PN expressions in fairly natural form, without worrying about calculating
all the variables needed for each expression, or manually accounting for the
various PN orders to which we may need to truncate PN expansions.  To achieve
this, we rely primarily on two macros.

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/macros.jl"]
Order   = [:macro, :module, :type, :constant, :function]
```


# Termination criteria

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


## Manipulating ODE solutions

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/combine_solutions.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


### Irrational constants

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


## Truncated series

We also have some utilities for dealing with series — or more
precisely summations, since we only handle finitely many terms.  In
particular, we need *truncated* multiplication (and some times
division) of truncated series.  This multiplication is associative and
there is a multiplicative identity element — though not always a
multiplicative inverse — which means that this structure naturally
forms a [monoid](https://en.wikipedia.org/wiki/Monoid).  (In fact, it
also forms more general structures, like a commutative algebra; all we
need is the monoidal structure.)

Here, we are assuming that there is a fixed order at which the series
are truncated, and that truncated multiplication preserves that order.
That is, if ``A`` and ``B`` are summations in terms up to ``v^N`` for
some integer ``N``, we want their product and ratio (if it exists) to
also be a summation in terms up to ``v^N``.

Note that we are not referring to these summations as "polynomials" in
``v``, because the coefficients will sometimes involve ``\ln(v)`` —
which is not technically permitted for polynomials.  In particular,
the presence of logarithms is irrelevant to our meaning of the "order"
of the truncated series.  This is standard practice in post-Newtonian
theory.[^1]

[^1]: Different texts in post-Newtonian theory treat these logarithmic
      terms with varying levels of rigor.  The preferred method is
      [Hadamard
      regularization](https://en.wikipedia.org/wiki/Hadamard_regularization)
      (often referred to in the literature as the *partie finie*).  A
      good summary is found in Section 6 of [Blanchet's Living
      Review](https://link.springer.com/article/10.12942/lrr-2014-2).
      Another potential approach could be taken following [this
      paper](https://doi.org/10.1016/0001-8708(89)90079-0).  But for
      our purposes, it will suffice to take the simplistic approach of
      treating logarithmic terms as if they were any other constant.

The following functions implement this behavior:

```@docs
PostNewtonian.truncated_series_inverse
PostNewtonian.truncated_series_product
PostNewtonian.truncated_series_ratio
```


## Miscellaneous

```@autodocs
Modules = [PostNewtonian]
Pages   = ["utilities/misc.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
