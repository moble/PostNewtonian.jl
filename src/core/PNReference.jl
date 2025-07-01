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
needs to be defined just in the given module.  One prominent example is that some authors
use `η` where we use `ν`.  In any case, this will make it easier for `@pn_expression` to
look through the current module for symbols that appear in the expression that need to be
replaced with a call to that function.

"""

function pn_reference(expr)
    # A module `Expr` has head `:module` and three arguments:
    # 1. a boolean indicating if this is a module (true) as opposed to a baremodule (false)
    # 2. the module name (a `Symbol`)
    # 3. the module body (an `Expr` with head `:block`)
    if MacroTools.isexpr(expr, :module)
        if length(expr.args) ≠ 3
            error(
                "Found a module expression with $(length(expr.args)) arguments:\n" *
                "$(expr)",
            )
        end
        if ! MacroTools.isexpr(expr.args[1], Bool)
            error("Found a module expression with a non-boolean flag: $(expr.args[1])")
        end
        if ! MacroTools.isexpr(expr.args[2], Symbol, :symbol)
            error(
                "Found a module expression with a non-symbol name: $(expr.args[2])\n" *
                "$(expr)",
            )
        end
        if ! MacroTools.isexpr(expr.args[3], :block)
            error("Found a module expression with a non-block body: $(dump(expr.args[3]))")
        end
        # Now, we assemble the new module, mostly by prepending some imports to the
        # contents
        new_module = Expr(
            :module,
            false,  # We turn this into a `baremodule`
            expr.args[2],  # This is the name of the module
            Expr(  # This is the new module body
                :block,
                :(import Base),
                :(eval(x) = Core.eval(Mod, x)),
                :(include(p) = Base.include(Mod, p)),
                :(using PostNewtonian: @pn_expression, @pn_expansion),
                expr.args[3].args...,  # The original module body
            ),
        )
        # Finally, we return this as a `:toplevel` expression, to ensure that it is
        # not buried in some block that causes an error.
        return Expr(:toplevel, new_module)
    else
        error("@pn_reference can only be used with module expressions, not `$expr`")
    end
end

"""
    @pn_reference

Create a module for a specific Post-Newtonian reference.

Different sources in the literature may use different notations for the same variables, so
enclosing definitions related to a particular reference in a module allows us to use the
appropriate notation within that module while still interfacing with this package and its
conventions outside of the module.  For example, some authors use `η` to denote the
symmetric mass ratio whereas this package uses `ν`.  We can import this with `using
..PostNewtonian: ν as η` inside the module, and then use `η` in the expressions defined
inside that module, to resemble the reference.  But when any expression using `η` is called
from outside the module, the choice of variable used inside that expression will not affect
the result.

This macro transforms the module to which it is applied, so that unlike a normal module that
effectively includes `using Base` — which allows access to such basic operations as
addition, multiplication, and so on — it instead includes `import Base`.  It also includes
`using PostNewtonian: @pn_expression, @pn_expansion`, the first of which redefines the basic
operations — at least those used in post-Newtonian expressions — to account for the number
type relevant to the input `PNSystem`.  By not importing the basic operations from `Base`,
we ensure that only those redefined by `@pn_expression` are available in the module, which
ensures that they preserve the number type.  If any are missing, `@pn_expression` should be
extended.

"""
@public macro pn_reference(ex)
    esc(pn_reference(ex))
end

@testitem "@pn_reference" begin
    using MacroTools: MacroTools
    using PostNewtonian: @pn_reference

    input = @macroexpand @pn_reference module Einstein1918
    using PostNewtonian: G, c, M, χ⃗₁, χ⃗₂, v, pn_order
    using Quaternionic: absvec
    @pn_expression χ₁(pnsystem) = absvec(χ⃗₁)
    @pn_expression χ₂(pnsystem) = absvec(χ⃗₂)
    const x = 3
    module InnerMod
        const z = 4
    end
    using .InnerMod: z
    @pn_expression function y(pnsystem)
        G*M/c^3 * @pn_expansion(x + χ₁ + χ₂ + z*(v/c))
    end
    end

    output = quote
        baremodule Einstein1918
        import Base
        eval(x) = Core.eval(Mod, x)
        include(p) = Base.include(Mod, p)
        using PostNewtonian: @pn_expression, @pn_expansion

        using PostNewtonian: G, c, M, χ⃗₁, χ⃗₂, v, pn_order
        using Quaternionic: absvec
        @pn_expression χ₁(pnsystem) = absvec(χ⃗₁)
        @pn_expression χ₂(pnsystem) = absvec(χ⃗₂)
        const x = 3
        module InnerMod
            const z = 4
        end
        using .InnerMod: z
        @pn_expression function y(pnsystem)
            G*M/c^3 * @pn_expansion(x + χ₁ + χ₂ + z*(v/c))
        end
        end
    end

    @test MacroTools.striplines(input).args[1] == MacroTools.striplines(output).args[1]
end
