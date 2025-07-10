public PNBase

# Documented below so that we can list the functions in the module automatically.
baremodule PNBase

# We want to be able to use the Base operators like `+`, `-`, `*`, `/`, and `^`, but we're
# already redefining those operators here, so we need alternative characters for them.
# Julia defines operators that can be used with equivalent precedence and as infix operators
# in this file: https://github.com/JuliaLang/julia/blob/master/src/julia-parser.scm.
# Circled versions of each will work, except for `^`, for which we use an arrow.  They can
# be typed as `\oplus`, `\ominus`, `\circledast`, `\oslash`, and `\uparrow`.
using Base: Base, + as ⊕, - as ⊖, * as ⊛, / as ⊘, ^ as ↑

using PostNewtonian: constant_convert, PNSystem
using PostNewtonian.InlineExports: @export
using Base: @inline

@export @inline ln(pnsystem::PNSystem, x) = Base.log(constant_convert(pnsystem, x))

@export @inline √(pnsystem::PNSystem, x) = Base.sqrt(constant_convert(pnsystem, x))

@export @inline (+)(pnsystem::PNSystem, x, y) = ⊕(constant_convert(pnsystem, x), y)
@inline (+)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
@inline function (+)(pnsystem::PNSystem, w, x, y, z...)
    return ⊕(constant_convert(pnsystem, w), x, y, z...)
end
@inline (+)(x::AbstractFloat, y::AbstractFloat) = x ⊕ y
@inline (+)(x::AbstractFloat) = ⊕(x)

@export @inline (-)(pnsystem::PNSystem, x, y) = ⊖(constant_convert(pnsystem, x), y)
@inline (-)(x::AbstractFloat, y::AbstractFloat) = x ⊖ y
@inline (-)(pnsystem::PNSystem, x) = ⊖(constant_convert(pnsystem, x))
@inline (-)(x::AbstractFloat) = ⊖(x)

@export @inline (*)(pnsystem::PNSystem, x, y) = ⊛(constant_convert(pnsystem, x), y)
@inline (*)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
@inline function (*)(pnsystem::PNSystem, w, x, y, z...)
    return ⊛(constant_convert(pnsystem, w), x, y, z...)
end
@inline (*)(x::AbstractFloat, y::AbstractFloat) = x ⊛ y

@export @inline (/)(pnsystem::PNSystem, x, y) = ⊘(constant_convert(pnsystem, x), y)
@inline (/)(x::AbstractFloat, y::AbstractFloat) = x ⊘ y

@export @inline (^)(pnsystem::PNSystem, x, n) = ↑(constant_convert(pnsystem, x), n)
@inline (^)(x::AbstractFloat, n::Integer) = x ↑ n
@inline (^)(x::AbstractFloat, n::AbstractFloat) = x ↑ n

# We can actually use the original definition of `//` for most purposes; we define it
# both ways here for simplicity.
@export @inline (//)(pnsystem::PNSystem, x, y) = Base.://(x, y)
@inline (//)(x, y) = Base.://(x, y)

end  # baremodule PNBase

const pnbase_functions = filter(s -> isa(getfield(PNBase, s), Function), names(PNBase))

@doc """
    PNBase

This module provides many of the same essential functions available in Julia's `Base`
module, but restricted to operations we expect to be used in post-Newtonian contexts, so
that we can ensure that the resulting expressions are valid within that framework.  This is
intended to be used inside `@pn_expression` modules.

It is intentionally very restrictive, so that it only includes the basic arithmetic
operations that are used in post-Newtonian expressions, and even then only in modified forms
that take a `PNSystem` as the first argument.  The defined operations are

    $(pnbase_functions)

This module is not intended to be used directly, but is imported by the
[`@pn_expression`](@ref) macro, and its methods are called by [`@pn_expansion`](@ref) to
ensure that the arithmetic operations preserve the number type of the input `PNSystem`.
"""
PNBase

@testitem "PNBase" begin
    baremodule Mod
    using PostNewtonian.PNBase
    end

    using PostNewtonian: pnbase_functions

    const ln = log

    for NT ∈ (Float16, Float64, BigFloat)
        pnsystem = BHNS(randn(NT, 15))
        z = (1, 2, 3, 17, 31, ℯ, π, 3//13, 47//59)

        for f ∈ (:+, :-, :*, :/, :^, :ln, :√)
            @test f ∈ pnbase_functions
        end
        for f ∈ (:+, :-, :*, :/, :^)
            for x ∈ z, y ∈ z
                @test eval(:(Mod.$f))(pnsystem, x, y) isa NT
                @test eval(:(Mod.$f))(pnsystem, x, y) == eval(:($f))(NT(x), NT(y))
            end
        end
        for f ∈ (:ln, :√)
            for x ∈ z
                @test eval(:(Mod.$f))(pnsystem, x) isa NT
                @test eval(:(Mod.$f))(pnsystem, x) == eval(:($f))(NT(x))
            end
        end
    end
end
