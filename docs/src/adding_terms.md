# Adding new PN terms or expressions

!!! danger "Cite your sources!"
    Remember that it is *absolutely crucial* to record the source of any
    expressions in the relevant docstrings, including a link to the relevant
    paper.  If the terms or expressions in a given function come from more than
    one source, use comments inside the code itself to clarify exactly which
    parts come from which source.

The first step in actually adding new information is to decide where
it fits in the hierarchy described in the ["Code structure" page](@ref
Structure-of-this-package's-code).  Most likely, you will want to add
new terms to existing PN expressions, and possibly simple functions in
[`PostNewtonian.DerivedVariables`](@ref "Derived variables").

Existing code should be a good guide on how to do this.  However, it
may appear as if some magic is happening in the various PN
expressions, whereby variables like `M₁`, `χ⃗₂`, `ν`, `Λ`, etc., can
be used without being defined.  These are automatically computed by
way of the [`@pn_expression`](@ref PostNewtonian.@pn_expression)
macro.  Also note that to correctly truncate a PN expansion at the
desired order, the [`@pn_expansion`](@ref PostNewtonian.@pn_expansion)
macro must be used.

To make it easier to compare code to original sources, we also want to
keep the code looking as similar as possible to equations as given in
those sources, so try to keep variable names the same, and order
things in the same ways as in the sources.  (Don't worry about
rewriting expressions to optimize evaluation; Julia will probably do a
better job automatically anyway.)  There are, however, a few important
exceptions to this rule:

   1. It is crucial to explicitly include factors of ``1/c`` to the
      appropriate power *relative to the 0-pN order term* inside the
      `@pn_expansion` macro.  For example, the first couple terms in
      the binding energy expansion look like
      ```julia
      32/5G * ν^2 * v^10 / c^5 * @pn_expansion(1 + (v/c)^2 * (-1247//336 - 35ν/12))
      ```
      The `1/c^5` factor at the beginning is useful to remind us of
      the *absolute* order of this expression, but the `1/c^2` factor
      inside the `@pn_expansion` macro is crucial to ensure that the
      expression can be correctly truncated at the desired order.
      Thus, if you want to add a new 6-pN term to the binding energy,
      you could either manually code it inside the `@pn_expansion`
      expression, including its `1/c^12` factor, or you could create a
      separate function that should look like this:
      ```julia
      @pn_expression function binding_energy_6pn(pnsystem::PNSystem)
          32/5G * ν^2 * v^12 / c^5 * @pn_expansion((v/c)^12 * 17//4)
      end
      ```
      Note that we have kept the `1/c^5` factor out front, and the
      `1/c^12` factor inside the `@pn_expansion` macro.  This tells
      the compiler that this is a 6-pN term, so if the user requests a
      lower pN order, the term should be ignored.
   2. It should be possible to evaluate PN expressions using different
      precisions.  To ensure this, enter fractions as `Irrational`s —
      e.g., `9//5` instead of `9/5`.  The latter would be immediately
      converted to the inexact 64-bit float value `1.8`, which would
      poison other values.  Note that if you are multiplying by
      something else that already has general float type, as in
      `9π/5`, you don't need to use `Irrational`; in this case, `9π`
      is evaluated first and achieves full precision, and is then
      divided by the exact integer 5, so that it retains full
      precision.  Note that a helpful regex to search for this case is
      `(?<!(ζ|n|\^))[0-9]+/[0-9]`, which has relatively few false
      positives in this repo.
   3. A slight caveat to the above is that an expression like `3λ/2`
      could *still* be converted to `Float64` if `λ` is defined at
      compile time to be `0`.  For type-stability reasons, Julia will
      always treat `0/2` just like it would treat `7/2`, which is
      converted to `Float64`.  Thus, it is probably safest to write
      expressions like `3//2 * λ`.
   4. If you happen to use any other math functions, similarly ensure
      that their arguments are converted appropriately to retain
      precision.  For unary functions, this can be done automatically
      by including the function name in the `unary_funcs` list used by
      [`@pn_expression`](@ref PostNewtonian.@pn_expression).
