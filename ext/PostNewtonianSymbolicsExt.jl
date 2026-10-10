module PostNewtonianSymbolicsExt

# See ../../src/predefinitions_Symbolics.jl for a few predefinitions of things that really
# only exist here, but will be needed elsewhere.  The documentation evidently needs to
# occur there as well.

using MacroTools: MacroTools
using SymbolicUtils: SymbolicUtils
isdefined(Base, :get_extension) ? (using Symbolics: Symbolics) : (import ..Symbolics)

using PostNewtonian
import PostNewtonian:
    hold,
    unhold,
    SymbolicPNSystem,
    type_converter,
    fundamental_quaternionic_variables,
    derived_variables,
    causes_domain_error!,
    prepare_pn_order,
    order_index,
    𝓔′,
    apply_to_first_add!,
    flatten_add!,
    flatten_mul!,
    pn_expression,
    pn_expansion,
    @pn_expansion,
    symbols,
    BBH, BHNS, NSNS,
    M₁,
    M₂,
    χ⃗₁,
    χ⃗₂,
    v,
    Φ,
    Λ₁,
    Λ₂,
    R,
    M,
    μ,
    ν,
    δ,
    q,
    ℳ,
    X₁,
    X₂,
    ln,
    ln2,
    ln3,
    ln5,
    ζ3,
    γₑ,
    _efficient_vector
#apply_to_first_add!, flatten_add!, pn_expression,
using RuntimeGeneratedFunctions: init, @RuntimeGeneratedFunction

using SciMLBase: parameterless_type

init(@__MODULE__)

function _efficient_vector(N, ::Type{Symbolics.Num})
    return Symbolics.variables(string(gensym()), 1:N)
end

### Moved from src/utilities/macros.jl

hold(x) = x
Symbolics.@register_symbolic hold(x)
Symbolics.derivative(::typeof(hold), args::NTuple{1,Any}, ::Val{1}) = 1

function unhold(expr)
    MacroTools.postwalk(expr) do x
        m = MacroTools.trymatch(:(f_(i_)), x)
        m === nothing || m[:f] !== hold ? x : Symbol(m[:i])
    end
end

function type_converter(::PNSystem{T}, x) where {T<:Symbolics.Num}
    return Symbolics.Num(SymbolicUtils.Term(hold, [x]))
end
function type_converter(::PNSystem{T}, x::Symbolics.Num) where {T<:Symbolics.Num}
    return x
end

# Add symbolic capabilities to all derived variables (fundamental variables already work)
for method ∈ [fundamental_quaternionic_variables; derived_variables]
    name = method.name
    @eval begin
        function PostNewtonian.$name(v::PNSystem{T}) where {T<:Symbolics.Num}
            return Symbolics.Num(SymbolicUtils.Sym{Real}(Symbol($name)))
        end
        function PostNewtonian.$name(v::Vector{T}) where {T<:Symbolics.Num}
            return Symbolics.Num(SymbolicUtils.Sym{Real}(Symbol($name)))
        end
    end
end

"""
    extract_var_factor(term, var)

Extract a factor of `var` from the product `term`.

This is a helper function for [`var_collect`](@ref).
"""
function extract_var_factor(term, var)
    if MacroTools.isexpr(term, :call) && term.args[1] ∈ ((/), :/)
        k₂, term₂ = extract_var_factor(term.args[2], var)
        k₃, term₃ = extract_var_factor(term.args[3], var)
        return k₂ - k₃, Expr(:call, term.args[1], term₂, term₃)
        #return k₂-k₃, :($(term.args[1]), $term₂, $term₃)
    end
    if !MacroTools.isexpr(term, :call) || term.args[1] ∉ ((*), :*)
        if term == var
            return 1, 1
        end
        m = MacroTools.trymatch(:((^)(v_, k_)), term)
        m = !isnothing(m) ? m : MacroTools.trymatch(:($(^)(v_, k_)), term)
        if !isnothing(m) && m[:v] == var
            return m[:k], 1
        end
        return 0, term
    end
    term = flatten_mul!(deepcopy(term))
    k = 0
    indices = Int[]
    for (i, factor) ∈ enumerate(term.args)
        if i == 1
            continue  # Skip the :*
        end
        if MacroTools.isexpr(factor, :call)
            k′, term′ = extract_var_factor(factor, var)
            # if term′ isa Expr
            #     term′ = Expr(:call, term′.args...)
            # end
            if k′ > 0
                k += k′
                term.args[i] = term′
            end
        else
            if factor == var
                k += 1
                push!(indices, i)
                continue
            end
            m = MacroTools.trymatch(:((^)(v_, k_)), factor)
            m = !isnothing(m) ? m : MacroTools.trymatch(:($(^)(v_, k_)), factor)
            if !isnothing(m) && m[:v] == var
                k += m[:k]
                push!(indices, i)
            end
        end
    end
    if !isempty(indices)
        splice!(term.args, indices)
    end
    # TODO: If there's just 1 arg left over, return it alone
    if length(term.args) == 2
        k, term.args[2]
    else
        k, term
    end
end

"""
    var_collect(expr, var)

Collect coefficients in `expr` of various powers of the variable `var`.

The inputs should be an
[`Expr`](https://docs.julialang.org/en/v1/manual/metaprogramming/#Program-representation)
and a single `Symbol` to be found in that `Expr`.

The return value is an `Int` representing the highest power of `v` in the expression, and an
`Expr` representing a `Tuple` of values corresponding to the coefficients of `var` to
various powers.  For example,
```jl-doctest
julia> PostNewtonian.var_collect(:(1 + a*v + b*v^2 + c*v^4), :v)
4, :((1, a, b, 0, c))
```
(Note that there was *no* factor in `v^3`.)  This result is convenient for passing to
`evalpoly`, for example.
"""
function var_collect(expr, var)
    if !MacroTools.isexpr(expr, :call)
        error("Input expression is not a call at its highest level: $expr")
    end
    terms = Dict{Int,Any}()
    if expr.args[1] ∉ ((+), :+, (-), :-)
        k, term = extract_var_factor(expr, var)
        terms[k] = term
    else
        for (i, term) ∈ enumerate(expr.args[2:end])
            k, term = extract_var_factor(term, var)
            if expr.args[1] ∈ ((-), :-) && i == 2
                if k ∈ keys(terms)
                    terms[k] = :($(terms[k]) - $(term))
                else
                    terms[k] = :(-$(term))
                end
            else
                if k ∈ keys(terms)
                    terms[k] = :($(terms[k]) + $(term))
                else
                    terms[k] = term
                end
            end
        end
    end
    max_k = maximum(keys(terms))
    term_exprs = [get(terms, k, 0) for k ∈ 0:max_k]
    return max_k, :(($(term_exprs...),))
end

function var_collect(expr::Symbolics.Num, var; max_power=100, max_gap=4)
    expr = SymbolicUtils.expand(expr)
    dict = Dict(var^j => 0 for j ∈ 1:max_power)
    c = SymbolicUtils.substitute(expr, dict; fold=false)
    expr = expr - c
    coefficients = [c]
    gap = 0
    for i ∈ 1:max_power
        dict[var ^ i] = 1
        if i > 1
            dict[var ^ (i - 1)] = 0
        end
        push!(coefficients, Symbolics.substitute(expr, dict; fold=false))
        if iszero(coefficients[end])
            gap += 1
            if gap ≥ max_gap
                return coefficients[1:(end - gap)]
            end
        else
            gap = 0
        end
    end
    return coefficients
end

macro pn_expansion(offset, pnsystem, expr)
    return esc(pn_expansion(offset, pnsystem, expr))
end

function pn_expansion(offset, pnsystem, expr)
    max_k, coefficients = var_collect(expr, :v)
    max_index = max_k + 1
    if offset != 0
        max_index_var = gensym("max_index")
        quote
            $max_index_var = min($max_index, order_index($pnsystem) - $offset)
            if $max_index_var < 1
                zero(v)
            else
                evalpoly(v, $coefficients[1:($max_index_var)])
            end
        end
    else
        :(evalpoly(v, $coefficients[1:min($max_index, order_index($pnsystem))]))
    end
end

## Moved from src/pn_systems.jl

causes_domain_error!(u̇, ::PNSystem{NT}) where {NT<:Symbolics.Num} = false

function SymbolicPNSystem(t::Type{PN}, PNOrder=typemax(Int)) where {PN <: QuasisphericalSystem}
    symbolic_vars = map(Symbolics.variable, symbols(PN))
    return PN{Symbolics.Num, Vector{Symbolics.Num}, prepare_pn_order(PNOrder)}(
        [symbolic_vars...]
    )
end

"""
    symbolic_pnsystem

TODO UPDATE

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

## Moved from src/pn_expressions/binding_energy.jl with a little padding to distinguish it
## from the new FastDifferentiation-based version

function make_symbolic_𝓔′(t::Type{PN}) where {PN <: QuasisphericalSystem}
    symbolic_pn = SymbolicPNSystem(t)
    let 𝓔 = 𝓔(symbolic_pn), v = v(symbolic_pn)
        ∂ᵥ = Symbolics.Differential(v)
        # Evaluate derivative symbolically
        𝓔′ = SymbolicUtils.simplify(Symbolics.expand_derivatives(∂ᵥ(𝓔)); expand=true)#, simplify_fractions=false)
        # Turn it into (an Expr of) a function taking one argument: `pnsystem`
        𝓔′ = Symbolics.build_function(𝓔′, :pnsystem; nanmath=false)
        # Remove `hold` (which we needed for Symbolics.jl to not collapse to Float64)
        𝓔′ = unhold(𝓔′)
        # "Flatten" the main sum, because Symbolics nests sums for some reason
        𝓔′ = apply_to_first_add!(𝓔′, flatten_add!)
        # Apply `@pn_expansion` to the main sum
        splitfunc = MacroTools.splitdef(𝓔′)
        splitfunc[:body] = apply_to_first_add!(splitfunc[:body], x -> :(@pn_expansion(-1, $x)))
        𝓔′ = MacroTools.combinedef(splitfunc)
        # Finally, apply the "macro" to it and get a full function out
        @RuntimeGeneratedFunction(pn_expression(1, 𝓔′))
    end
end



const 𝓔′Symbolics =
    Dict(pn => make_symbolic_𝓔′(pn)
         for pn ∈ (BBH, BHNS, NSNS))

function 𝓔′(
    pnsystem::PNSystem{NT,ST,PNOrder},
    ::Val{:Symbolics};
    pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum),
) where {NT,ST,PNOrder,PNExpansionReducer}
    if PNExpansionReducer != sum
        error("Symbolic 𝓔′ is not implemented for PNExpansionReducer other than `sum`")
    else
        𝓔′Symbolics[parameterless_type(pnsystem)](pnsystem)
    end
end

end #module
