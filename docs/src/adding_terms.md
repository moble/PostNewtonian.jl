# Adding new PN terms or expressions

The first step in actually adding new information is to decide where it fits in
the hierarchy above.  Most likely, you will want to add new terms to existing PN
expressions (item 4 above).  Existing code should be a good guide on how to do
this, but note that if you will be using derived variables inside your
expressions that don't yet exist inside this package, you should define
functions that take exactly one `PNSystem` argument in the
`DerivedVariables` submodule.

Remember that it is *absolutely crucial* to record the source of any expressions
in the relevant docstrings, including a link to the relevant paper.  If the
terms or expressions in a given function come from more than one source, use
comments inside the code itself to clarify exactly which parts come from which
source.

We also want to keep the code looking as similar as possible to equations given
in original sources, so try to keep variable names the same, and order things in
the same ways as in the sources.  (Don't worry about rewriting expressions to
optimize evaluation; Julia will probably do a better job automatically anyway.)
There are, however, a few important exceptions to this rule:

   1. We use a PN-expansion parameter `v`, instead of `x`.  Consistency with
      this is important to simplify our treatment of PN expressions
      symbolically.  This isn't too difficult, because `v ↔ √x` is easy to do
      visually.  Just remember that `ln(x) = 2ln(v)`.
   2. Due to problems in defining a Taylor series in `v` when the coefficients
      include factors like `ln(v)`, we have to treat `ln(v)` as a constant
      symbolically (even though it will get updated automatically at each step).
      So just use `lnv` instead.
   3. It should be possible to evaluate PN expressions using different
      precisions.  To ensure this, enter fractions as `Irrational`s — e.g.,
      `3//2` instead of `3/2`.  The latter would be immediately converted to the
      64-bit float `1.5`, which would poison other values.  Note that if you are
      multiplying by something else that already has general float type, as in
      `3π/2`, you don't need to use `Irrational`.  In this case, `3π` is
      evaluated first and achieves full precision, and is then divided by the
      exact integer 2, so that it retains full precision.
   4. If you happen to use any other math functions, similarly ensure that their
      arguments are converted appropriately to retain precision.  For unary
      functions, this can be done automatically by including the function name
      in the `unary_funcs` list used by
      [`@compute_pn_variables`](@ref PostNewtonian.@compute_pn_variables).