# Structure of this package's code

## Code hierarchy

There is a fairly simple hierarchy to the code:

1. **System and state**

   Objects of type `AbstractPNSystem` represent the PN system itself and some information
   about it, including:
   - The `Binary` type — `BBH`, `BHNS`, or `NSNS`
   - The float type `T` — `Float64`, etc.
   - The PN expansion order `PNOrder` — a `Rational`
   - The `Expansion` type — `TaylorT1`, `TaylorT4`, or `TaylorT5`
   - The current `state` of the system, including the fundamental variables

   The first of these is represented as the type itself, the last is stored as a
   vector inside the object, and the rest are represented as type parameters.

   Note that basically everything below will be written as a function of such an
   object, which we will denote `pnsystem`.

2. **Fundamental variables** 
   
   This consists of the basic variables describing the system such as `M₁`,
   `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`.  For systems with matter, this may also include
   Love numbers `λ₁` and `λ₂`.

   It's important to note that these should all be accessed through functions
   like `M₁(pnsystem)` rather than directly like `pnsystem.M₁`.  This allows
   Julia's type system to get involved, enabling important optimizations.

   Also, these variables can be automatically computed in functions that need
   them with the `@compute_pn_variables` macro.  For example, you can directly
   use the symbols `M₁`, `M₂`, etc., in a function that is wrapped in that
   macro, without any qualifiers to specify where those variables are coming
   from, and the macro will automatically and efficiently evaluate them for you
   because they are defined in the `PostNewtonian.FundamentalVariables` module.

3. **Derived variables**

   These are variables that are frequently used in PN expressions that can be
   computed from the fundamental variables, and are *defined* in terms of them.
   For example, the total mass `M ≔ M₁+M₂`, the antisymmetric spin vector `χ⃗ₐ ≔
   (χ⃗₁-χ⃗₂)/2`, or the orbital separation vector `n̂ ≔ R x̂ R̄`.  While none of
   these are strictly necessary, it is helpful to be able to write the same
   variables used in articles providing the PN expressions in the code itself.

   Because they are defined solely in terms of fundamental variables, which can
   be computed from an `AbstractAbstractPNSystem` alone, these are all written
   as functions of such an object — such as `M(pnsystem)`, `χ⃗ₐ(pnsystem)`, and
   `n̂(pnsystem)`.

   Again, these quantities will be automatically computed for you in any
   function wrapped in the `@compute_pn_variables` macro because they are
   defined in the `PostNewtonian.DerivedVariables` module.

4. **PN expressions**

   Unlike derived variables, these are not *defined* in terms of the fundamental
   variables, but they can be calculated in terms of both fundamental and
   derived variables.  These are generally the result of post-Newtonian
   expansions — the most important examples being the flux [`𝓕`](@ref), binding
   energy [`𝓔`](@ref), and the waveform mode weights [`h!`](@ref) themselves.

5. **Evaluation**

   This is where the ODE integration actually occurs, to evolve the orbital
   dynamics of the system, and the computation of the waveform mode weights in
   terms of those dynamics.

6. **Compatibility layers**

   This is an optional level of abstraction that allows us to wrap the
   evaluation layer in functions that are designed to look and act more like
   other packages' functions.  As of this writing, the only such layer is for
   [`GWFrames`](https://github.com/moble/GWFrames) compatibility, but similar
   wrappers could certainly be added.

## Adding new PN terms or expressions

The first step in actually adding new information is to decide where it fits in
the hierarchy above.  Most likely, you will want to add new terms to existing PN
expressions (item 4 above).  Existing code should be a good guide on how to do
this, but note that if you will be using derived variables inside your
expressions that don't yet exist inside this package, you should define
functions that take exactly one `AbstractPNSystem` argument in the
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
      in the `unary_funcs` list used by [`@compute_pn_variables`](@ref).
