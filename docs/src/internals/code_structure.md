# Structure of this package's code

There is a fairly simple hierarchy to the code.  Beyond some basic utilities,
which are related to how we write the code, rather than post-Newtonian
calculations *per se*, we have the following:

1. **System and state**

   Objects of type `PNSystem` represent the PN system itself and some
   information about it, including:
   - The `Binary` type — `BBH`, `BHNS`, or `NSNS`
   - The float type `T` — `Float64`, etc.
   - The PN expansion order `PNOrder` — a `Rational`
   - The current `state` of the system, including the fundamental variables

   The first of these is represented as the type itself, the last is stored as a
   vector inside the object, and the rest are represented as type parameters.

   Note that basically everything below will be written as a function of such an
   object, which we will denote `pnsystem`.

2. **Fundamental variables** 
   
   This consists of the basic variables describing the system such as `M₁`,
   `M₂`, `χ⃗₁`, `χ⃗₂`, `R`, `v`.  For systems with matter, this may also include
   tidal deformability for each star, `Λ₁` and `Λ₂`.

   It's important to note that these should all be accessed through functions
   like `M₁(pnsystem)` rather than directly like `pnsystem.M₁`.  This allows
   Julia's type system to get involved, enabling important optimizations.

   Also, these variables can be automatically computed in functions that need
   them with the `@pn_expression` macro.  For example, you can directly
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
   be computed from an `PNSystem` alone, these are all written
   as functions of such an object — such as `M(pnsystem)`, `χ⃗ₐ(pnsystem)`, and
   `n̂(pnsystem)`.

   Again, these quantities will be automatically computed for you in any
   function wrapped in the `@pn_expression` macro because they are
   defined in the `PostNewtonian.DerivedVariables` module.

4. **PN expressions**

   Unlike derived variables, these are not *defined* in terms of the fundamental
   variables, but they can be calculated in terms of both fundamental and
   derived variables.  These are generally the result of post-Newtonian
   expansions — the most important examples being the flux [`𝓕`](@ref), binding
   energy [`𝓔`](@ref), and the waveform mode weights [`h!`](@ref) themselves.

5. **PN expansions**

   These are the specific parts of a "PN expression" that are expanded
   in powers of ``1/c``.  For example, the binding energy `𝓔` is
   expanded in powers of `v/c`.  When PN expansions are defined using
   the [`@pn_expression`](@ref PostNewtonian.@pn_expression) macro,
   the expansion is automatically truncated at the appropriate order.
   See [`PNTerm`](@ref PostNewtonian.PNTerm) and [`PNExpansion`](@ref
   PostNewtonian.PNExpansion) for more details.

6. **Dynamics**

   This is where the ODE integration actually occurs, to evolve the orbital
   dynamics of the system.

7. **Evaluation**

   Finally, we construct the waveforms themselves.  This level contains the main
   interface that will usually be used from Julia code, and should be restricted
   to fairly high-level functions like `PNWaveform(M₁, M₂, ...)`, while still
   handling the full range of options that will be present in "Dynamics", for
   example.

8. **Compatibility layers**

   This is an optional level of abstraction that allows us to wrap the
   evaluation layer in functions that are designed to look and act more like
   other packages' functions.  As of this writing, the only such layer is for
   [`GWFrames`](https://github.com/moble/GWFrames) compatibility, but similar
   wrappers could certainly be added.
