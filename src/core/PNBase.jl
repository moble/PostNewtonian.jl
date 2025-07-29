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

using PostNewtonian: constant_convert, PNSystem, ExactNumber, ExactIntegerBased
using PostNewtonian.InlineExports: @export
using Base: @inline

@export @inline ln(pnsystem::PNSystem, x) = Base.log(constant_convert(pnsystem, x))
@inline ln(x::AbstractFloat) = Base.log(x)

@export @inline √(pnsystem::PNSystem, x) = Base.sqrt(constant_convert(pnsystem, x))
@inline √(x::AbstractFloat) = Base.sqrt(x)

@export @inline (+)(pnsystem::PNSystem, x, y) = ⊕(constant_convert(pnsystem, x), y)
@inline (+)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
# @inline function (+)(pnsystem::PNSystem, w, x, y, z...)
#     return ⊕(constant_convert(pnsystem, w), x, y, z...)
# end
@inline (+)(x::AbstractFloat, y::AbstractFloat) = x ⊕ y
@inline (+)(x::AbstractFloat, y::ExactNumber) = x ⊕ y
@inline (+)(x::ExactNumber, y::AbstractFloat) = x ⊕ y
@inline (+)(x::ExactIntegerBased, y::ExactIntegerBased) = x ⊕ y
@inline (+)(x::Number) = x
(+)(x, y, z, args...) = Base.afoldl(+, (x + y) + z, args...)

@export @inline (-)(pnsystem::PNSystem, x, y) = ⊖(constant_convert(pnsystem, x), y)
@inline (-)(x::AbstractFloat, y::AbstractFloat) = x ⊖ y
@inline (-)(x::AbstractFloat, y::ExactNumber) = x ⊖ y
@inline (-)(x::ExactNumber, y::AbstractFloat) = x ⊖ y
@inline (-)(x::ExactIntegerBased, y::ExactIntegerBased) = x ⊖ y
@inline (-)(x::Number) = ⊖(x)
@inline (-)(pnsystem::PNSystem, x) = ⊖(constant_convert(pnsystem, x))

@export @inline (*)(pnsystem::PNSystem, x, y) = ⊛(constant_convert(pnsystem, x), y)
@inline (*)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
# @inline function (*)(pnsystem::PNSystem, w, x, y, z...)
#     return ⊛(constant_convert(pnsystem, w), x, y, z...)
# end
@inline (*)(x::AbstractFloat, y::AbstractFloat) = x ⊛ y
@inline (*)(x::AbstractFloat, y::ExactNumber) = x ⊛ y
@inline (*)(x::ExactNumber, y::AbstractFloat) = x ⊛ y
@inline (*)(x::ExactIntegerBased, y::ExactIntegerBased) = x ⊛ y
(*)(x, y, z, args...) = Base.afoldl(*, (x * y) * z, args...)

@export @inline (/)(pnsystem::PNSystem, x, y) = ⊘(constant_convert(pnsystem, x), y)
@inline (/)(x::AbstractFloat, y::AbstractFloat) = x ⊘ y
@inline (/)(x::AbstractFloat, y::ExactNumber) = x ⊘ y
@inline (/)(x::ExactNumber, y::AbstractFloat) = x ⊘ y

@export @inline (^)(pnsystem::PNSystem, x, n) = ↑(constant_convert(pnsystem, x), n)
@inline (^)(x::AbstractFloat, n::ExactNumber) = x ↑ n
@inline (^)(x::ExactNumber, n::AbstractFloat) = x ↑ n
@inline (^)(x::AbstractFloat, n::AbstractFloat) = x ↑ n

# We can actually use the original definition of `//` for most purposes; we define it
# both ways here for simplicity.
@export @inline (//)(pnsystem::PNSystem, x, y) = Base.://(x, y)
@inline (//)(x, y) = Base.://(x, y)

end  # baremodule PNBase

# COV_EXCL_START  (Covered, but missed due to const-folding)
const pnbase_functions = filter(s -> isa(getfield(PNBase, s), Function), names(PNBase))
# COV_EXCL_END

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
        Z = (1, 2, 3, 17, 31, 3//13, 47//59)
        z = (Z..., ℯ, π)

        for x ∈ z
            @test eval(:(Mod.:+))(pnsystem, x) isa NT
            @test eval(:(Mod.:+))(pnsystem, x) == NT(x)
            @test eval(:(Mod.:+))(x) == x

            @test eval(:(Mod.:-))(pnsystem, x) isa NT
            @test eval(:(Mod.:-))(pnsystem, x) == -NT(x)
            @test eval(:(Mod.:-))(x) == -x

            @test eval(:(Mod.:*))(pnsystem, x) isa NT
            @test eval(:(Mod.:*))(pnsystem, x) == NT(x)
        end

        for f ∈ (:+, :*)
            op = f≡:+ ? sum : prod
            for w ∈ Z, x ∈ Z, y ∈ Z
                @test eval(:(Mod.$f))(pnsystem, w, x, y) isa NT
                @test eval(:(Mod.$f))(pnsystem, w, x, y) == op(NT, (w, x, y))
                @test eval(:(Mod.$f))(w, x, y) ≡ op((w, x, y))
                for v ∈ Z
                    @test eval(:(Mod.$f))(pnsystem, v, w, x, y) isa NT
                    @test eval(:(Mod.$f))(pnsystem, v, w, x, y) == op(NT, (v, w, x, y))
                    @test eval(:(Mod.$f))(v, w, x, y) ≡ op((v, w, x, y))
                end
            end
        end

        for f ∈ (:+, :-, :*, :/, :^, :ln, :√)
            @test f ∈ pnbase_functions
        end
        for f ∈ (:+, :-, :*, :/, :^)
            for x ∈ z, y ∈ z
                @test eval(:(Mod.$f))(pnsystem, x, y) isa NT
                @test eval(:(Mod.$f))(pnsystem, x, y) == eval(:($f))(NT(x), NT(y))
                @test eval(:(Mod.$f))(x, NT(y)) isa NT
                @test eval(:(Mod.$f))(x, NT(y)) ≈ eval(:($f))(NT(x), NT(y))
                @test eval(:(Mod.$f))(NT(x), y) isa NT
                @test eval(:(Mod.$f))(NT(x), y) ≈ eval(:($f))(NT(x), NT(y))
                @test eval(:(Mod.$f))(NT(x), NT(y)) isa NT
                @test eval(:(Mod.$f))(NT(x), NT(y)) == eval(:($f))(NT(x), NT(y))
            end
        end
        for f ∈ (:+, :-, :*, ://)
            for x ∈ Z, y ∈ Z
                @test typeof(eval(:(Mod.$f))(x, y)) ≡ typeof(eval(:($f))(x, y))
                @test eval(:(Mod.$f))(x, y) == eval(:($f))(x, y)
            end
        end
        for x ∈ Z, y ∈ Z
            @test typeof(eval(:(Mod.://))(pnsystem, x, y)) ≡ typeof(x // y)
            @test eval(:(Mod.://))(pnsystem, x, y) == x // y
        end
        for f ∈ (:ln, :√)
            for x ∈ z
                @test eval(:(Mod.$f))(pnsystem, x) isa NT
                @test eval(:(Mod.$f))(pnsystem, x) == eval(:($f))(NT(x))
                @test eval(:(Mod.$f))(NT(x)) == eval(:($f))(NT(x))
            end
            for x ∈ Z
                @test_throws MethodError eval(:(Mod.$f))((x))
            end
        end
    end
end
