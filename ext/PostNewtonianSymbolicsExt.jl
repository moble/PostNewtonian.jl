module PostNewtonianSymbolicsExt

# See ../../src/predefinitions_Symbolics.jl for a few predefinitions of things that really
# only exist here, but will be needed elsewhere.  The documentation evidently needs to
# occur there as well.

import MacroTools
import SymbolicUtils
isdefined(Base, :get_extension) ? (import Symbolics) : (import ..Symbolics)

using PostNewtonian
import PostNewtonian: hold, unhold, SymbolicPNSystem,
    type_converter, fundamental_quaternionic_variables, derived_variables,
    causes_domain_error!, prepare_pn_order,
    apply_to_first_add!, flatten_add!, pn_expression, order_index,
    M₁, M₂, χ⃗₁, χ⃗₂, v, Φ, Λ₁, Λ₂,
    R, M, μ, ν, δ, q, ℳ, X₁, X₂,
    ln, ln2, ln3, ln5, ζ3, γₑ,
    _efficient_vector


function _efficient_vector(N, ::Type{Symbolics.Num})
    Symbolics.variables(string(gensym()), 1:N)
end

### Moved from src/utilities/macros.jl

hold(x) = x
Symbolics.@register_symbolic hold(x)
Symbolics.derivative(::typeof(hold), args::NTuple{1,Any}, ::Val{1}) = 1

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


function SymbolicPNSystem(PNOrder=typemax(Int))
    Symbolics.@variables M₁ M₂ χ⃗₁ˣ χ⃗₁ʸ χ⃗₁ᶻ χ⃗₂ˣ χ⃗₂ʸ χ⃗₂ᶻ Rʷ Rˣ Rʸ Rᶻ v Φ Λ₁ Λ₂
    ET = typeof(M₁)
    SymbolicPNSystem{Vector{ET}, prepare_pn_order(PNOrder), ET}(
        [M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v, Φ],
        Λ₁, Λ₂
    )
end


"""
    symbolic_pnsystem

A symbolic `PNSystem` that contains symbolic information for all types of `PNSystem`s.

In particular, note that this object has an (essentially) infinite `PNOrder` and has nonzero
values for quantities like `Λ₁` and `Λ₂`.  If you want different choices, you may need to
call [`SymbolicPNSystem`](@ref) yourself, or even construct a different specialized subtype
of `PNSystem` (it's not hard).

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


# ## Moved from src/pn_expressions/binding_energy.jl and renamed
# const 𝓔′Symbolics = let 𝓔=𝓔(symbolic_pnsystem), v=v(symbolic_pnsystem)
#     ∂ᵥ = Symbolics.Differential(v)
#     # Evaluate derivative symbolically
#     𝓔′ = SymbolicUtils.simplify(Symbolics.expand_derivatives(∂ᵥ(𝓔)), expand=true)#, simplify_fractions=false)
#     # Turn it into (an Expr of) a function taking one argument: `pnsystem`
#     𝓔′ = Symbolics.build_function(𝓔′, :pnsystem, nanmath=false)
#     # Remove `hold` (which we needed for Symbolics.jl to not collapse to Float64)
#     𝓔′ = unhold(𝓔′)
#     # "Flatten" the main sum, because Symbolics nests sums for some reason
#     𝓔′ = apply_to_first_add!(𝓔′, flatten_add!)
#     # Apply `@pn_expansion` to the main sum
#     splitfunc = MacroTools.splitdef(𝓔′)
#     splitfunc[:body] = apply_to_first_add!(
#         splitfunc[:body],
#         x->:(@pn_expansion(-1, $x))
#     )
#     𝓔′ = MacroTools.combinedef(splitfunc)
#     # Finally, apply the "macro" to it and get a full function out
#     eval(pn_expression(1, 𝓔′))::Function
# end

end #module
