# Pre-define a few functions / structs that the Symbolics extension can extend

"""
    hold(x)

Delay evaluation of the argument in `Symbolics` expressions.

This is just a helper function that acts trivially — like the `identity` function — but also
gets registered with `Symbolics` to avoid evaluation of the argument.  For example, we can
preserve expressions like `π^2`, which Julia would normally convert directly to a `Float64`.

Note that you probably don't want to use this function directly; this will probably be done
for you by [`@pn_expression`](@ref PostNewtonian.@pn_expression) or similar.  If you *do*
want to use this directly, you probably want another layer of indirection to construct
something like `Symbolics.Num(SymbolicUtils.Term(hold, [x]))` so that you can use the result
in a symbolic expression.
"""
function hold end

"""
    unhold(expr)

Remove occurrences of [`hold`](@ref) from an `Expr`.
"""
function unhold end

"""
    SymbolicPNSystem{ST, PNOrder, ET}(state, Λ₁, Λ₂)

A `PNSystem` that contains information as variables from
[`Symbolics.jl`](https://symbolics.juliasymbolics.org/).

# Examples
```jldoctest
julia> using Symbolics

julia> using PostNewtonian: M₁, M₂, χ⃗₁, χ⃗₂, SymbolicPNSystem

julia> symbolic_pnsystem = SymbolicPNSystem()
SymbolicPNSystem{Vector{Num}, 9223372036854775805//2, Num}(Num[M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ], Λ₁, Λ₂)

julia> M₁(symbolic_pnsystem), M₂(symbolic_pnsystem)
(M₁, M₂)

julia> χ⃗₁(symbolic_pnsystem)
χ⃗₁

julia> χ⃗₂(symbolic_pnsystem)
χ⃗₂
```
"""
struct SymbolicPNSystem{NT,ST,PNOrder} <: PNSystem{NT,ST,PNOrder}
    state::ST
    Λ₁::NT
    Λ₂::NT
end
