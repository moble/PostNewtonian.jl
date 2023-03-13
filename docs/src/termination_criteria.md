# Termination criteria for inspiral evolution

Hopefully, it should not be necessary to directly use these termination
criteria.  They are used by default in the [`inspiral`](@ref) function.  But
certain particularly extreme physical parameters may lead ODEs that are
difficult to integrate — especially if new PN systems or terms are introduced.
These, or similar functions may be helpful examples of
["callbacks"](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
that can be passed to the ODE integrator.

```@meta
CurrentModule = PostNewtonian
```

```@docs
termination_forwards
termination_backwards
dtmin_terminator
nonfinite_terminator
```
