# Adding new PN terms or expressions

!!! danger "Cite your sources!"
    Remember that it is *absolutely crucial* to record the source of any
    expressions in the relevant docstrings, including a link to the relevant
    paper.  If the terms or expressions in a given function come from more than
    one source, use comments inside the code itself to clarify exactly which
    parts come from which source.

The first step in actually adding new information is to decide where it fits in
the hierarchy described in the ["Code structure" page](@ref
Structure-of-this-package's-code).  Most likely, you will want to add new terms
to existing PN expressions, and possibly simple functions in
[`PostNewtonian.DerivedVariables`](@ref "Derived variables").

Existing code should be a good guide on how to do this.  However, it may appear
as if some magic is happening in the various PN expressions, whereby variables
like `M₁`, `χ⃗₂`, `ν`, `Λ`, etc., can be used without being defined.  These are
automatically computed by way of the [`@pn_expression`](@ref
PostNewtonian.@pn_expression) macro.  Also note that to correctly truncate a PN
expansion at the desired order, the [`@pn_expansion`](@ref
PostNewtonian.@pn_expansion) macro must be used.

To make it easier to compare code to original sources, we also want to keep the
code looking as similar as possible to equations as given in those sources, so
try to keep variable names the same, and order things in the same ways as in the
sources.  (Don't worry about rewriting expressions to optimize evaluation; Julia
will probably do a better job automatically anyway.)  There are, however, a few
important exceptions to this rule:

   1. We use a PN-expansion parameter `v`, instead of `x`.  Consistency with
      this is important to simplify our treatment of PN expressions
      symbolically.  This isn't too difficult, because `v ↔ √x` is easy to do
      visually.  Just remember that `ln(x) = 2ln(v)`.
   2. It should be possible to evaluate PN expressions using different
      precisions.  To ensure this, enter fractions as `Irrational`s — e.g.,
      `3//2` instead of `3/2`.  The latter would be immediately converted to the
      64-bit float `1.5`, which would poison other values.  Note that if you are
      multiplying by something else that already has general float type, as in
      `3π/2`, you don't need to use `Irrational`; in this case, `3π` is
      evaluated first and achieves full precision, and is then divided by the
      exact integer 2, so that it retains full precision.
   3. A slight caveat to the above is that an expression like `3λ/2` could
      *still* be converted to `Float64` if `λ` is defined at compile time to be
      `0`.  For type-stability reasons, Julia will always treat `0/2` just like
      it would treat `7/2`, which is converted to `Float64`.  Thus, it is
      probably safest to write expressions like `3//2 * λ`.
   4. If you happen to use any other math functions, similarly ensure that their
      arguments are converted appropriately to retain precision.  For unary
      functions, this can be done automatically by including the function name
      in the `unary_funcs` list used by
      [`@pn_expression`](@ref PostNewtonian.@pn_expression).
