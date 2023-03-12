# Structure of this package's code

## Code hierarchy

There is a fairly simple hierarchy to the code:

1. **State**

   Objects of type `PNState` represent the PN system itself and some information
   about it, including:
   - The float type `T` — `Float64`, etc.
   - The PN expansion order `PNOrder` — a `Rational`
   - The `Binary` type — `BBH`, `BHNS`, or `NSNS`
   - The `Expansion` type — `TaylorT1`, `TaylorT4`, or `TaylorT5`
   - The current state of the system, including the fundamental variables

   All but the last item are stored as type parameters; the last is stored as a vector inside the object.

   Note that basically everything below will be written as a function of such an
   object, which we will denote `pnstate`.

2. **Fundamental variables** 
   
   This consists of the basic variables describing the system such as `M₁`,
   `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`.  For systems with matter, this may also include
   Love numbers `λ₁` and `λ₂`.

   It's important to note that these should all be accessed through functions
   like `M₁(pnstate)` rather than directly like `pnstate.M₁`.  This allows
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
   be computed from a `PNState` alone, these are all written as functions of
   such an object — such as `M(pnstate)`, `χ⃗ₐ(pnstate)`, and `n̂(pnstate)`.

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

## Adding new PN terms or expressions

The first step in actually adding new information is to decide where it fits in
the hierarchy above.  Most likely, you will want to add new terms to existing PN
expressions (item 4 above).  Existing code should be a good guide on how to do
this, but note that if you will be using derived variables inside your
expressions that don't yet exist inside this package, you should define
functions that take exactly one `PNState` argument in the `DerivedVariables`
submodule.

Remember that it is *absolutely crucial* to record the source of any expressions
in the relevant docstrings, including a link to the relevant paper.  If the
terms or expressions in a given function come from more than one source, use
comments inside the code itself to clarify exactly which parts come from which
source.
