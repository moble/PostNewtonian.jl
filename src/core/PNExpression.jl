macro pn_expression(arg_index, func=:(nothing))
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
    # default value `Val(sum)`.  The `Val` is captured in the "where" parameter
    # `PNExpansionReducer`
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

    # # Next, we walk the expression tree and find any function calls to one of the
    # # `pnexpressionarithmetic_functions`, and insert `pnsystem` as the first argument.
    new_body = MacroTools.postwalk(new_body) do x
        if iscall(x, pnexpressionarithmetic_functions)
            insert!(x.args, 2, (pnsystem))
            x
        else
            x
        end
    end

    # Finally, include the `let` statements so that the functions can all be used as plain
    # symbols without having to call them on `pnsystem`, but do so without burying it in a
    # code block.
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
        f * χ₁ˣ
    end
    @pn_expression function j(pnsyst)
        f / χ₁ʸ
    end
    @pn_expression function k(pnsyste)
        f ^ χ₁ᶻ
    end
    @pn_expression function l(pnsystm)
        f * ln(v)
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
    end

    # Test that `f` and only `f` is the thing visible to the macro function finder
    # Test the output against BBH{Float16}, BBH{Float64}, and BBH{BigFloat}

    # @pn_expression function test_arithmetic(pnsystem)
    #     return (pnsystem.M + pnsystem.c) * (pnsystem.M - pnsystem.c) / pnsystem.c^2
    # end
    # result = test_arithmetic(pnsystem)
    # @test result == (pnsystem.M^2 - pnsystem.c^2) / pnsystem.c^2

end
