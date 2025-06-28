## Possibilities for Code Structure

### `StaticArrays.FieldVector`

This is a nice idea.  Given a `struct` with a series of fields, making
it a subtype of `FieldVector` makes it an `SVector` or `MVector`
(depending on whether it is immutable or mutable).  The fields can be
accessed by the names you give them, but everything else is like a
nicely completed `StaticArray`, which is handled well by all the SciML
packages.  So, for example, we would define `BBH` as

```julia
using StaticArrays: FieldVector

struct BBH{T, PNOrder} <: FieldVector{14, T}
    M₁::T
    M₂::T
    χ⃗₁::T
    χ⃗₁ˣ::T
    χ⃗₁ʸ::T
    χ⃗₁ᶻ::T
    χ⃗₂::T
    χ⃗₂ˣ::T
    χ⃗₂ʸ::T
    χ⃗₂ᶻ::T
    R::T
    Rʷ::T
    Rˣ::T
    Rʸ::T
    Rᶻ::T
    v::T
    Φ::T
end
```

The big disadvantage is that it is not possible to contain non-state
variables in the `struct` — for example, the tidal coupling
parameters.  Making them into state variables would make the ODE
system larger than it needs to be.  Making them into type parameters
would be bad design, lengthening compile times.  I don't see any other
way to incorporate them into the `struct`, so this option is
discarded.

### `LabelledArrays`

This is a nice design, and works well with SciML.  But it doesn't
allow any other metadata — even `PNOrder`, which would be possible
with `FieldVector`.  So this is discarded.

### Emulating `LabelledArrays`' use of `Syms`

Part of the type of a `LabelledArray` is a `Syms` object, which is a
`Tuple` of `Symbol`s.  I could basically reimplement this for
`PNSystem`s, which could bring a great deal of flexibility.  For
example, if a `PNSystem` has an `e` field, it must be an eccentric
system.  So, in principle, we could dispatch on `Syms`.  This seems
clunky to me.  Also, it would remove our ability to further
specialize, unless we defined adequate abstract super types.  I don't
have a solid argument against this, but I think it would be better to
just define a `symbols` function for each subtype.

### Make `symbols` as function of type, and emulate `LabelledArrays` features

This will be a lot of work, but basically copying all the methods of
`LabelledArrays` seems like the way to go.

### Separate types for mutable and immutable systems

I think it would be too much work to maintain two separate types for
each of the `PNSystem`s, one mutable and one immutable.  Probably
better to just store the container type as a type parameter, and test
whether that container `ismutabletype` — for example, inside
`setindex!`, we could test and raise a more informative error if the
system `!ismutabletype`.

## Diagram

```mermaid
flowchart TB
  %% define each layer as its own box
  subgraph Core["<code>core</code><br/>Building blocks of the code"]
    <!-- a["XYZ"]
    b["CYS"] -->
  end
  subgraph Systems["<code>pn_systems</code><br/>Types encoding various binaries"]
  end
  subgraph Literature["<code>literature</code><br/>Modules for each reference containing pieces of PN expressions"]
  end
  subgraph Expressions["<code>pn_expressions</code><br/>Functions for computing physical quantities"]
  end
  subgraph Interface["<code>interface</code><br/>High-level functions for users to call"]
  end

  %% draw the arrows in build‐order
  Core --> Systems
  Systems --> Literature
  Literature --> Expressions
  Expressions --> Interface
```

- PostNewtonian
  - interface: high-level functions for users to call
    - orbital_evolution
    - waveform
    - pn
  - pn_expressions: functions for computing physical quantities
    - flux
    - energy
    - angular_momentum
    - precession
    - dynamics
    - waveforms
  - literature: modules for each reference containing pieces of PN expressions
    - common variables
    - Ref1
      - import fundamental variables
      - import common variables
      - define variables not in common variables
      - write individual expressions as separate functions with `@pn_expression`
  - pn_systems: types encoding various binaries
    - PNSystem
      - AbstractBBHSystem
        - QuasicircularBBH
        - QuasisphericalBBH
        - EccentricNonspinningBBH
        - BBH
      - AbstractBHNSSystem
        - BHNS
      - AbstractNSNSSystem
        - NSNS
  - core: the building blocks of the code
    - PNExpression
    - PNReference
    - PNExpansion
    - PNTerm

- `PNSystem`
  - Define abstract subtypes:
    - `AbstractBBHSystem`
      - `QuasicircularBBH`
      - `QuasisphericalBBH`
      - `EccentricNonspinningBBH`

- `PNExpression`
- `PNExpansion`
- `PNTerm`

- literature
  - common variables
  - each reference
    - directory and module named by bibtex key
    - modules are imported directly under top-level `PostNewtonian`
    - explicitly import any fundamental variables used in the expressions
    - import variables that are defined the same as the common variables
    - define variables that are not defined in the common variables
    - write individuals

- PN expressions
  - Assemble functions for high-level quantities like flux, etc.
  - Inside each function, simply call the expressions from the
    literature, prepending with module name.

`@pn_expression` will

1. look for symbols in the expression, and for any that matches a function
   name in the current module, add a `let` binding for that symbol to equal
   the function called on the `pnsystem` argument.
