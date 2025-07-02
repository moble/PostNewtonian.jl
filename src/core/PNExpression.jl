"""
    @pn_expression [arg_index=1] func

This macro makes it easier to define post-Newtonian expressions.

The optional first argument to this macro is `arg_index`, which just tells us which argument
to the function `func` is a `PNSystem` — defaulting to the first argument.  For example, the
variables defined in [`PostNewtonian.FundamentalVariables`](@ref "Fundamental variables")
all take a single argument of `pnsystem`, which is used to compute the values for those
variables; this macro just needs to know where to find `pnsystem`.

The required argument to this macro is a function the will compute some values related to
the `pnsystem`.

Once it has this information, there are three types of transformations the macro will make:

 1. Adds a keyword argument `pn_expansion_reducer::Val{PNExpansionReducer}=Val(sum)` to the
    function signature.  This is used to determine how to reduce the PN expansion terms.
    The default is `Val(sum)`, which will just return a single number,  but `Val(identity)`
    can be used to return the expansion.  The reducing function can be used inside the main
    function as `PNExpansionReducer`, and will be automatically used inside any
    `@pn_expansion`.
 2. For every function defined in the module where this macro is called that takes a single
    `PNSystem` argument, this macro will look for the name of that function inside the
    expression.  If it finds any such names, it will replace them with the corresponding
    calls to the function with `pnsystem` as the argument.  For example, you can simply use
    the symbols `M₁` or `μ` in your code, rather than calling them as [`M₁(pnsystem)`](@ref
    M₁) or [`μ(pnsystem)`](@ref μ) every time they appear.  To be more explicit, this is
    achieved by defining the relevant quantities in a `let` block placed around the body of
    `func`, so that the values may be used efficiently without recomputation.  If you need
    to use one of those functions with different arguments, you'll need to address them
    explicitly with the module name — as in `PostNewtonian.v(;Ω, M)`.
 3. Insert the `pnsystem` argument as the first argument to each occurrence of
    `@pn_expansion` that needs it.

For example, we might write a function like this:
```julia
@pn_expression function f(pnsystem)
    5μ/c^2 * @pn_expansion(1 + 4(v/c)^2)
end
```
That is effectively rewritten by this macro as
```julia
function f(pnsystem, ::Val{PNExpansionReducer}=Val(sum)) where {PNExpansionReducer}
    let μ=μ(pnsystem), c=c(pnsystem), v=v(pnsystem)
        5μ/c^2 * @pn_expansion pnsystem (1 + 4(v/c)^2)
    end
end
```
The [`@pn_expansion`](@ref) macro is further expanded so that the final result looks like
```julia
function f(pnsystem, ::Val{PNExpansionReducer}=Val(sum)) where {PNExpansionReducer}
    let μ=μ(pnsystem), c=c(pnsystem), v=v(pnsystem)
        5μ/c^2 * (
            let c=PNExpansionParameter(pnsystem)
                PNExpansionReducer(1 + 4(v/c)^2)
            end
        )
    end
end
```
Moreover, this expression should itself be defined within an [`@pn_reference`](@ref) module,
which preserves the precision of the `pnsystem`, so that we can define `pnsystem` with,
e.g., `Float16` or `BigFloat` numbers, and natural expressions like `2π` or `4/3` won't
automatically be converted to `Float64` and spoil the requested precision.
"""
@public macro pn_expression(arg_index, func=:(nothing))
    # Note that macros secretly get two extra variables: `__module__` and `__source__`.
    # These refer to the place where the macro was called.  The second is useful for
    # debugging, but the first is crucial here, because it lets us see what's defined in the
    # module in which this macro was used.

    # First, we just search for everything defined in the module where this macro is called.
    all_names = if VERSION < v"1.12.0-beta1"
        @warn "Some `names` in `@pn_expression` at $(__source__)\n" *
            "may not be found in Julia 1.11 or earlier; for now, be sure to use\n" *
            "`import` rather than `using` inside `@pn_reference` modules." maxlog=10
        names(__module__; all=true, imported=true)
    else
        names(__module__; all=true, imported=true, usings=true)
    end

    # Next, we narrow that list down to callable objects defined in the `__module__` where
    # this macro was used.  Specifically, we filter for those objects that take one
    # `PNSystem`, and nothing else, as an argument.  These will be valid targets for
    # replacement in the body that `@pn_expression` is applied to.
    pnsystem_functions = filter(
        function (name)
            if isdefined(__module__, name)
                f = getfield(__module__, name)
                if isa(f, Base.Callable)
                    if !isempty(methods(f, Tuple{PNSystem}))
                        return true
                    end
                end
            end
            return false
        end, all_names
    )

    # Now, we pass that information along to the function that actually modifies our code
    return esc(pn_expression(arg_index, func, pnsystem_functions, __module__, __source__))
end

function pn_expression(arg_index, func, pnsystem_functions, __module__, __source__)
    # Handle the default argument
    if func === :(nothing)
        arg_index, func = 1, arg_index
    end

    # Check to make sure that `arg_index` is an integer; we'll check that it's valid below.
    if !isa(arg_index, Integer)
        error(
            "The optional argument to `@pn_expression` must be an integer, not " *
            "$(typeof(arg_index))\nThe incorrect application is at $(__source__).",
        )
    end

    # Check to make sure that `func` is a function-definition expression
    if !MacroTools.isdef(func)
        error(
            "The `@pn_expression` macro can only be applied to function definitions.\n" *
            "The incorrect application is at $(__source__).",
        )
    end

    # MacroTools has this nice function to parse function definitions so we can operate on
    # each piece separately.
    splitfunc = MacroTools.splitdef(func)

    # Check that the argument index is valid
    if arg_index < 1 || arg_index > length(splitfunc[:args])
        error(
            "Invalid argument index $arg_index for function $(splitfunc[:name]); " *
            "it must be between 1 and $(length(splitfunc[:args]))\n" *
            "The incorrect application is at $(__source__).",
        )
    end

    # Get the name of the PNSystem argument
    pnsystem = MacroTools.namify(splitfunc[:args][arg_index])

    # Append a keyword argument `pn_expansion_reducer` to the function signature, with
    # default value `Val(sum)`.  The `Val` parameter is captured as `PNExpansionReducer`,
    # which we can then use in the body of the function — specifically, it is automatically
    # used by the `@pn_expansion` macro.
    splitfunc[:kwargs] = [
        splitfunc[:kwargs]
        :($(Expr(:kw, :(pn_expansion_reducer::Val{PNExpansionReducer}), :(Val(Base.sum)))))
    ]

    # Add `PNExpansionReducer` to the `where` clause (see previous comment)
    splitfunc[:whereparams] = (splitfunc[:whereparams]..., :PNExpansionReducer)

    # Now look for any of the `pnsystem_functions` in the body of the function and prepare
    # an entry for the `let` block that will surround the body on output.
    pnsystem_function_exprs = [
        :($s = $s($pnsystem)) for
        s ∈ filter(s -> MacroTools.inexpr(splitfunc[:body], s), pnsystem_functions)
    ]

    # Next, add `pnsystem` as the argument to each @pn_expansion call...
    new_body = MacroTools.postwalk(splitfunc[:body]) do x
        if MacroTools.isexpr(x, :macrocall) &&
            x.args[1] == Symbol("@pn_expansion") &&
            !isa(x.args[end - 1], Symbol)  # ...unless there's already one there
            x′ = deepcopy(x)
            insert!(x′.args, length(x′.args), pnsystem)
            x′
        else
            x
        end
    end

    # Now expand all the macros inside the body
    new_body = macroexpand(__module__, new_body; recursive=true)

    # Next, we walk the expression tree and find all function calls to any of the
    # `pnexpressionarithmetic_functions`, and insert `pnsystem` as the first argument.
    new_body = MacroTools.postwalk(new_body) do x
        if iscall(x, pnexpressionarithmetic_functions)
            insert!(x.args, 2, (pnsystem))
            x
        else
            x
        end
    end

    # Finally, include the `let` statements so that the functions can all be used as plain
    # symbols without having to call them on `pnsystem`, but we do so without burying it in
    # a code block.
    new_body = MacroTools.unblock(quote
        #@fastmath
        let $(pnsystem_function_exprs...)
            $(MacroTools.unblock(new_body))
        end
    end)

    # Replace the old body with the new one
    splitfunc[:body] = new_body

    # And reassemble the function definition
    return MacroTools.combinedef(splitfunc)
end

@testitem "@pn_expression" begin
    baremodule Mod
    # These are the expressions added by `@pn_reference`
    using Base: Base, Val
    eval(x::Expr) = Core.eval(Base.@__MODULE__, x)
    include(p::AbstractString) = Base.include(Base.@__MODULE__, p)
    using PostNewtonian: PostNewtonian, @pn_expression, @pn_expansion, 𝒾, γₑ, ζ3
    using PostNewtonian.PNExpressionArithmetic

    # Below are what we would normally write manually
    import PostNewtonian: M₁, M₂, χ₁ˣ, χ₁ʸ, χ₁ᶻ, v, Φ
    const x = Base.://(7, 2)
    f(pnsys::PostNewtonian.PNSystem) = x
    @pn_expression function g(pnsy)
        f + M₁
    end
    @pn_expression h(pns) = f - M₂
    @pn_expression function i(pn)
        f * χ₁ˣ * γₑ
    end
    @pn_expression function j(pnsyst)
        f / χ₁ʸ
    end
    @pn_expression function k(pnsyste)
        f ^ χ₁ᶻ
    end
    @pn_expression function l(pnsystm)
        f * ln(v) * ζ3
    end
    @pn_expression function m(pnsystem)
        √f * Φ
    end
    end  # baremodule Mod

    import .Mod

    const ln = log

    for NT ∈ (Float16, Float64, BigFloat)
        pnsystem = BHNS(NT.(1:15))
        for v ∈ (:g, :h, :i, :j, :k, :l, :m)
            @test eval(:(Mod.$v))(pnsystem) isa NT
        end
        @test Mod.g(pnsystem) == NT(7//2) + pnsystem[:M₁]
        @test Mod.h(pnsystem) == NT(7//2) - pnsystem[:M₂]
        @test Mod.i(pnsystem) == NT(7//2) * pnsystem[:χ₁ˣ] * PostNewtonian.γₑ
        @test Mod.j(pnsystem) == NT(7//2) / pnsystem[:χ₁ʸ]
        @test Mod.k(pnsystem) == NT(7//2) ^ pnsystem[:χ₁ᶻ]
        @test Mod.l(pnsystem) == NT(7//2) * log(pnsystem[:v]) * PostNewtonian.ζ3
        @test Mod.m(pnsystem) == √(NT(7//2)) * pnsystem[:Φ]
    end
end
