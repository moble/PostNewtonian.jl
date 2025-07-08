public PNBase

"""
    PNBase

This module provides many of the same essential functions available in Julia's `Base`
module, but restricted to operations we expect to be used in post-Newtonian contexts, so
that we can ensure that the resulting expressions are valid within that framework.
"""
baremodule PNBase
import Base
using PostNewtonian: constant_convert, PNSystem, PNTerm, PNExpansionParameter, PNExpansion
using PostNewtonian.InlineExports: @export
using Base: @inline

@export @inline ln(pnsystem::PNSystem, x) = Base.log(constant_convert(pnsystem, x))

@export @inline √(pnsystem::PNSystem, x) = Base.sqrt(constant_convert(pnsystem, x))

@export @inline (+)(pnsystem::PNSystem, x, y) = Base.:+(constant_convert(pnsystem, x), y)
@inline (+)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
@inline function (+)(pnsystem::PNSystem, w, x, y, z...)
    return Base.:+(constant_convert(pnsystem, w), x, y, z...)
end

@export @inline (-)(pnsystem::PNSystem, x, y) = Base.:-(constant_convert(pnsystem, x), y)

@export @inline (*)(pnsystem::PNSystem, x, y) = Base.:*(constant_convert(pnsystem, x), y)
@inline (*)(pnsystem::PNSystem, x) = constant_convert(pnsystem, x)
@inline function (*)(pnsystem::PNSystem, w, x, y, z...)
    return Base.:*(constant_convert(pnsystem, w), x, y, z...)
end

@export @inline (/)(pnsystem::PNSystem, x, y) = Base.:/(constant_convert(pnsystem, x), y)

@export @inline (^)(pnsystem::PNSystem, x, n) = Base.:^(constant_convert(pnsystem, x), n)

end  # baremodule PNBase

const pnbase_functions = filter(s -> isa(getfield(PNBase, s), Function), names(PNBase))

@doc """
    PNBase

This module provides arithmetic operations to be used inside `@pn_expression` modules.

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
