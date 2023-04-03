# Derived variables

## Orbital elements

```@docs
n̂
λ̂
ℓ̂
Ω
```

## Mass combinations

```@docs
PostNewtonian.M
PostNewtonian.μ
PostNewtonian.ν
PostNewtonian.δ
PostNewtonian.q
PostNewtonian.ℳ
```

## Spin combinations

```@docs
S⃗₁
S⃗₂
S⃗
Σ⃗
χ⃗
χ⃗ₛ
χ⃗ₐ
χₑ
χₚ
```

Additionally, there are numerous convenience functions to give certain
components of the spins.  They take a single `pnsystem` argument and are not
exported.  Given the definitions above, they are all fairly self explanatory —
such as `χ₁²`, which gives `χ⃗₁ ⋅ χ⃗₁`; or `χ₁₂ = χ⃗₁ ⋅ χ⃗₂`; or `Sₙ = S⃗ ⋅ n̂`.
Like all the other fundamental and derived variables, these can be used directly
in PN expressions modified by the [`@pn_expression`](@ref
PostNewtonian.@pn_expression) macro.
