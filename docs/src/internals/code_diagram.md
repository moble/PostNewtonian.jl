# Diagram

```mermaid
flowchart TB
  %% define each layer as its own box
  subgraph Core["<code>core</code><br/>Building blocks of the code"]
    a["XYZ"]
    b["CYS"]
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
