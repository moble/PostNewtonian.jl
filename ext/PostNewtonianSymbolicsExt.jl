module PostNewtonianSymbolicsExt

using PostNewtonian
import PostNewtonian: type_converter, fundamental_quaternionic_variables, derived_variables,
    var_collect, causes_domain_error!, prepare_pn_order,
    apply_to_first_add!, flatten_add!, pn_expression,
    M, μ, ν, δ, q, ℳ, X₁, X₂,
    M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ, Λ₁, Λ₂

using MacroTools
using SymbolicUtils
isdefined(Base, :get_extension) ? (using Symbolics) : (using ..Symbolics)

export SymbolicPNSystem, symbolic_pnsystem, var_collect

### Moved from src/utilities/macros.jl

"""
    hold(x)

Delay evaluation of the argument in `Symbolics` expressions.

This is just a helper function that acts trivially — like the `identity` function — but also
gets registered with `Symbolics` to avoid evaluation of the argument.  For example, we can
preserve expressions like `π^2`, which Julia would normally convert directly to a `Float64`.

Note that you probably don't want to use this function directly; this will probably be done
for you by [`@pn_expression`](@ref) or similar.  If you *do* want to use this directly, you
probably want another layer of indirection to construct something like
`Symbolics.Num(SymbolicUtils.Term(hold, [x]))` so that you can use the result in a symbolic
expression.
"""
hold(x) = x
Symbolics.@register_symbolic hold(x)
Symbolics.derivative(::typeof(hold), args::NTuple{1,Any}, ::Val{1}) = 1

"""
    unhold(expr)

Remove occurrences of [`hold`](@ref) from an `Expr`.
"""
function unhold(expr)
    MacroTools.postwalk(expr) do x
        m = MacroTools.trymatch(:(f_(i_)), x)
        m === nothing || m[:f]!==hold ? x : Symbol(m[:i])
    end
end

function type_converter(::PNSystem{T}, x) where {T<:Vector{Symbolics.Num}}
    Symbolics.Num(SymbolicUtils.Term(hold, [x]))
end
function type_converter(::PNSystem{T}, x::Symbolics.Num) where {T<:Vector{Symbolics.Num}}
    x
end

# Add symbolic capabilities to all derived variables (fundamental variables already work)
for method ∈ [fundamental_quaternionic_variables; derived_variables]
    name = method.name
    @eval begin
        function PostNewtonian.$name(v::PNSystem{T}) where {T<:Vector{Symbolics.Num}}
            Symbolics.Num(SymbolicUtils.Sym{Real}(Symbol($name)))
        end
        function PostNewtonian.$name(v::Vector{T}) where {T<:Symbolics.Num}
            Symbolics.Num(SymbolicUtils.Sym{Real}(Symbol($name)))
        end
    end
end

function var_collect(expr::Symbolics.Num, var; max_power=100, max_gap=4)
    expr = SymbolicUtils.expand(expr)
    dict = Dict(var^j => 0 for j=1:max_power)
    c = SymbolicUtils.substitute(expr, dict, fold=false)
    expr = expr - c
    coefficients = [c]
    gap = 0
    for i in 1:max_power
        dict[var^i] = 1
        if i > 1
            dict[var^(i-1)] = 0
        end
        push!(coefficients, Symbolics.substitute(expr, dict, fold=false))
        if iszero(coefficients[end])
            gap += 1
            if gap ≥ max_gap
                return coefficients[1:end-gap]
            end
        else
            gap = 0
        end
    end
    coefficients
end


## Moved from src/systems.jl

causes_domain_error!(u̇, ::PNSystem{VT}) where {VT<:Vector{Symbolics.Num}} = false


"""
    SymbolicPNSystem{T, PNOrder}(state, Λ₁, Λ₂)

A `PNSystem` that contains information as variables from
[`Symbolics.jl`](https://symbolics.juliasymbolics.org/).

See also [`symbolic_pnsystem`](@ref) for a particular general instance of this type.
"""
struct SymbolicPNSystem{ST, PNOrder, ET} <: PNSystem{ST, PNOrder}
    state::ST
    Λ₁::ET
    Λ₂::ET

    function SymbolicPNSystem(PNOrder=typemax(Int))
        Symbolics.@variables M₁ M₂ χ⃗₁ˣ χ⃗₁ʸ χ⃗₁ᶻ χ⃗₂ˣ χ⃗₂ʸ χ⃗₂ᶻ Rʷ Rˣ Rʸ Rᶻ v Φ Λ₁ Λ₂
        ET = typeof(M₁)
        new{Vector{ET}, prepare_pn_order(PNOrder), ET}(
            [M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ],
            Λ₁, Λ₂
        )
    end
end

"""
    symbolic_pnsystem

A symbolic `PNSystem` that contains symbolic information for all types of `PNSystem`s.

In particular, note that this object has an (essentially) infinite `PNOrder`, uses the
`TaylorT1` approximant, and has nonzero values for quantities like `Λ₁` and `Λ₂`.  If you
want different choices, you may need to call [`SymbolicPNSystem`](@ref) yourself, or even
construct a different specialized subtype of `PNSystem` (it's not hard).

# Examples
```jldoctest
julia> using PostNewtonian: M₁, M₂, χ⃗₁, χ⃗₂

julia> M₁(symbolic_pnsystem), M₂(symbolic_pnsystem)
(M₁, M₂)

julia> χ⃗₁(symbolic_pnsystem)
χ⃗₁

julia> χ⃗₂(symbolic_pnsystem)
χ⃗₂
```
"""
const symbolic_pnsystem = SymbolicPNSystem()


## Moved from src/fundamental_variables.jl
Λ₁(pn::SymbolicPNSystem) = pn.Λ₁
Λ₂(pn::SymbolicPNSystem) = pn.Λ₂


## Moved from src/pn_expressions/binding_energy.jl and renamed
const 𝓔′Symbolics = let 𝓔=𝓔(symbolic_pnsystem), v=v(symbolic_pnsystem)
    ∂ᵥ = Symbolics.Differential(v)
    # Evaluate derivative symbolically
    𝓔′ = SymbolicUtils.simplify(Symbolics.expand_derivatives(∂ᵥ(𝓔)), expand=true)#, simplify_fractions=false)
    # Turn it into (an Expr of) a function taking one argument: `pnsystem`
    𝓔′ = Symbolics.build_function(𝓔′, :pnsystem, nanmath=false)
    # Remove `hold` (which we needed for Symbolics.jl to not collapse to Float64)
    𝓔′ = unhold(𝓔′)
    # "Flatten" the main sum, because Symbolics nests sums for some reason
    𝓔′ = apply_to_first_add!(𝓔′, flatten_add!)
    # Apply `@pn_expansion` to the main sum
    splitfunc = MacroTools.splitdef(𝓔′)
    splitfunc[:body] = apply_to_first_add!(
        splitfunc[:body],
        x->:(@pn_expansion(-1, $x))
    )
    𝓔′ = MacroTools.combinedef(splitfunc)
    # Finally, apply the "macro" to it and get a full function out
    eval(pn_expression(1, 𝓔′))::Function
end

end #module
