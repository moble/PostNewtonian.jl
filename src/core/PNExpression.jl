public PNExpressionArithmetic

baremodule PNExpressionArithmetic
import Base
#using ..PostNewtonian: constant_convert, PNSystem
using PostNewtonian: constant_convert, PNSystem
using PostNewtonian.InlineExports: @export
using Base: @inline

@export @inline ln(pnsystem::PNSystem, x) = Base.log(constant_convert(pnsystem, x))
@export const log = ln

@export @inline √(pnsystem::PNSystem, x) = Base.sqrt(constant_convert(pnsystem, x))
@export const sqrt = √

@export @inline (+)(pnsystem::PNSystem, x, y) = Base.:+(constant_convert(pnsystem, x), y)
@export @inline (-)(pnsystem::PNSystem, x, y) = Base.:-(constant_convert(pnsystem, x), y)
@export @inline (*)(pnsystem::PNSystem, x, y) = Base.:*(constant_convert(pnsystem, x), y)
@export @inline (/)(pnsystem::PNSystem, x, y) = Base.:/(constant_convert(pnsystem, x), y)
@export @inline (^)(pnsystem::PNSystem, x, n) = Base.:^(constant_convert(pnsystem, x), n)

@export @inline function (+)(pnsystem::PNSystem, w, x, y, z...)
    return Base.:+(constant_convert(pnsystem, w), x, y, z...)
end
@export @inline function (*)(pnsystem::PNSystem, w, x, y, z...)
    return Base.:*(constant_convert(pnsystem, w), x, y, z...)
end

end  # baremodule PNExpressionArithmetic

const pnexpressionarithmetic_functions = filter(
    s -> isa(getfield(PNExpressionArithmetic, s), Function), names(PNExpressionArithmetic)
)

@doc """
    PNExpressionArithmetic

This module provides arithmetic operations to be used inside `@pn_expression` modules.

It is intentionally very restrictive, so that it only includes the basic arithmetic
operations that are used in post-Newtonian expressions, and even then only in modified forms
that take a `PNSystem` as the first argument.  The defined operations are

    $(pnexpressionarithmetic_functions)

This module is not intended to be used directly, but is imported by the
[`@pn_expression`](@ref) macro, and its methods are called by [`@pn_expansion`](@ref) to
ensure that the arithmetic operations preserve the number type of the input `PNSystem`.
"""
PNExpressionArithmetic

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
                    if !isempty(methods(f, Tuple{PNSystem}, __module__))
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
            @show x.head x.args
            println()
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
            $(new_body)
        end
    end)

    # Replace the old body with the new one
    splitfunc[:body] = new_body

    # And reassemble the function definition
    return MacroTools.combinedef(splitfunc)
end

# baremodule Mod
# import Base
# (+)(x, y) = Base.:+(Base.BigFloat(x), y)
# const x = Base.π + 3
# end
