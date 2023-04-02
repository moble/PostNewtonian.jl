For convenience, an assortment of constructors is provided for relevant classes
of physical scenarios, as well as the ability to construct [random](@ref rand)
systems across a broad range of reasonable physical parameters.

In each case, the returned object is a `PNSystem`, which can be used as input to
most functions in this package — most notably the `orbital_evolution` function.
For example, to integrate the inspiral of a binary black-hole system in
["superkick"](@ref superkick) configuration, we could write
```julia
inspiral = orbital_evolution(superkick())
```
This implicitly provides the `M₁`, `M₂`, `χ⃗₁`, `χ⃗₂`, `Ωᵢ`, `Rᵢ`, and `PNOrder`
arguments (along with `λ₁` and/or `λ₂` for `BHNS` or `NSNS` systems) to
[`orbital_evolution`](@ref); other keyword arguments to that function can be
provided after `superkick()`.



```@autodocs
Modules = [PostNewtonian]
Pages   = ["assorted_binaries/examples.jl"]
```

```@docs
rand(pnclass::PNSystem; v, PNOrder)
```