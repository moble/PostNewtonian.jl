"""
This file defines the `@pn_reference` macro, which is used to create a module for a specific
Post-Newtonian reference.  Its purpose is to *un*define even the most basic operations, such
as addition, multiplication, and so on, to ensure that all operations are sufficiently
overridden by `@pn_expression`.  For example, integer division `1/3` will result in a
`Float64` by default.  But `@pn_expression` will extract the number type `NT` of the
`pnsystem` argument, and redefine integer division to return something compatible with `NT`
instead — possibly a `Double64` or `Sympy` expression, for example.

The [julia docs
say](https://docs.julialang.org/en/v1/manual/modules/#Default-top-level-definitions-and-bare-modules)
that a standard module is essentially this:

```julia
baremodule Mod

using Base

eval(x) = Core.eval(Mod, x)
include(p) = Base.include(Mod, p)

...

end
```
We just want to replace that `using` with `import`.  This will ensure that only operations
that are already handled by `@pn_expression` will be defined; if they aren't, there will be
an error.

Similarly, we want the author to *have to* define every single variable to be used in the PN
expression.  In most cases, the variable can just be imported from `common_variables` under
its own name.  There may be cases, however, where the name needs to change, or a variable
needs to be defined just in the given module.  In any case, this will make it easier for
`@pn_expression` to look through the current module for symbols that appear in the
expression that need to be replaced with a call to that function.


"""
