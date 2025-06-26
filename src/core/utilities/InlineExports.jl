# NOTE: This file is a modified version of the `InlineExports.jl` package, which is
#       licensed under the MIT License. The original source can be found at:
#       https://github.com/dalum/InlineExports.jl

module InlineExports

import Base: @__doc__

export @export, @public

eval(quote
    """
        @export

    Return the expression with all bindings exported.

    ```
    julia> module M
               using InlineExports
               @export begin
                   const a = 2
                   abstract type S <: Number end
                   struct T <: S
                       val
                   end
               end
               @export f(x::TT) where {TT<:S} = x.val^2
           end
    M

    julia> using .M

    julia> f(T(a))
    4
    ```
    """
    macro $(Symbol("export"))(expr::Expr)
        r = handle(expr)
        if r isa Symbol
            return quote
                export $(esc(r))
                @__doc__ $(esc(expr))
            end
        else
            return quote
                export $(map(esc, r)...)
                @__doc__ $(esc(expr))
            end
        end
    end
end)

eval(quote
    """
        @public

    Return the expression with all bindings marked as public.

    ```
    julia> module M
               using InlineExports
               @public begin
                   const a = 2
                   abstract type S <: Number end
                   struct T <: S
                       val
                   end
               end
               @public f(x::TT) where {TT<:S} = x.val^2
           end
    M

    julia> using .M

    julia> f(T(a))
    4
    ```
    """
    macro $(Symbol("public"))(expr::Expr)
        r = handle(expr)
        if r isa Symbol
            return quote
                public $ (esc(r))
                @__doc__ $(esc(expr))
            end
        else
            return quote
                public $ (map(esc, r)...)
                @__doc__ $(esc(expr))
            end
        end
    end
end)

handle(::Any) = nothing
handle(x::Symbol) = x
handle(expr::Expr) = handle(Val(expr.head), expr)

handle(::Val{:block}, expr) = filter(x -> x !== nothing, map(handle, expr.args))
handle(::Val{:const}, expr) = handle(expr.args[1])
handle(::Val{:(::)}, expr) = handle(expr.args[1])
handle(::Val{:(=)}, expr) = handle(expr.args[1])
handle(::Val{:function}, expr) = handle(expr.args[1])
handle(::Val{:where}, expr) = handle(expr.args[1])
handle(::Val{:macro}, expr) = Symbol("@", handle(expr.args[1]))
handle(::Val{:struct}, expr) = handle(expr.args[2])
handle(::Union{Val{:abstract},Val{:primitive}}, expr) = handle(expr.args[1])

handle(::Val{:<:}, expr) = handle(expr.args[1])
handle(::Val{:curly}, expr) = handle(expr.args[1])
handle(::Val{:call}, expr) = handle(expr.args[1])
handle(::Val{:macrocall}, expr) = filter(x -> x !== nothing, map(handle, expr.args[3:end]))

end  # module InlineExports

@testitem "InlineExports" begin
    module Bla
    using ..InlineExports: @export, @public

    @export begin
        "`const a` doc"
        const a = 2
        "`abstract type S` doc"
        abstract type S <: Number end
        "`struct T` doc"
        struct T <: S
            val
        end
    end
    """
        f(x)

    Here a doc!
    """
    @export f(x::TT) where {TT<:S} = x.val^2

    @public begin
        "`const b` doc"
        const b = 3
        "`abstract type U` doc"
        abstract type U <: Number end
        "`struct V` doc"
        struct V <: U
            val
        end
    end
    """
        g(x)

    There a doc!
    """
    @public g(x::VV) where {VV<:U} = x.val^3

    begin
        "`const c` doc"
        const c = 5
        "`abstract type W` doc"
        abstract type W <: Number end
        "`struct X` doc"
        struct X <: W
            val
        end
    end
    """
        h(x)

    Everywhere a doc, doc!
    """
    h(x::XX) where {XX<:W} = x.val^5

    end

    using .Bla

    @test f(T(a)) == 4
    @test Bla.g(Bla.V(Bla.b)) == 27
    @test Bla.h(Bla.X(Bla.c)) == 125

    @test Base.isexported(Bla, :a)
    @test Base.isexported(Bla, :S)
    @test Base.isexported(Bla, :T)
    @test Base.isexported(Bla, :f)
    @test !Base.isexported(Bla, :b)
    @test !Base.isexported(Bla, :U)
    @test !Base.isexported(Bla, :V)
    @test !Base.isexported(Bla, :g)
    @test !Base.isexported(Bla, :c)
    @test !Base.isexported(Bla, :W)
    @test !Base.isexported(Bla, :X)
    @test !Base.isexported(Bla, :h)

    @test Base.ispublic(Bla, :a)
    @test Base.ispublic(Bla, :S)
    @test Base.ispublic(Bla, :T)
    @test Base.ispublic(Bla, :f)
    @test Base.ispublic(Bla, :b)
    @test Base.ispublic(Bla, :U)
    @test Base.ispublic(Bla, :V)
    @test Base.ispublic(Bla, :g)
    @test !Base.ispublic(Bla, :c)
    @test !Base.ispublic(Bla, :W)
    @test !Base.ispublic(Bla, :X)
    @test !Base.ispublic(Bla, :h)

    @test Base.hasproperty(Bla, :a)
    @test Base.hasproperty(Bla, :S)
    @test Base.hasproperty(Bla, :T)
    @test Base.hasproperty(Bla, :f)
    @test Base.hasproperty(Bla, :b)
    @test Base.hasproperty(Bla, :U)
    @test Base.hasproperty(Bla, :V)
    @test Base.hasproperty(Bla, :g)
    @test Base.hasproperty(Bla, :c)
    @test Base.hasproperty(Bla, :W)
    @test Base.hasproperty(Bla, :X)
    @test Base.hasproperty(Bla, :h)

    @test Base.@__doc__ Bla.a == "`const a` doc"
    @test Base.@__doc__ Bla.b == "`const b` doc"
    @test Base.@__doc__ Bla.c == "`const c` doc"
    @test Base.@__doc__ Bla.S == "`abstract type S` doc"
    @test Base.@__doc__ Bla.U == "`abstract type U` doc"
    @test Base.@__doc__ Bla.W == "`abstract type W` doc"
    @test Base.@__doc__ Bla.T == "`struct T` doc"
    @test Base.@__doc__ Bla.V == "`struct V` doc"
    @test Base.@__doc__ Bla.X == "`struct X` doc"
    @test Base.@__doc__ Bla.f == "f(x)\n\nHere a doc!"
    @test Base.@__doc__ Bla.g == "g(x)\n\nThere a doc!"
    @test Base.@__doc__ Bla.h == "h(x)\n\nEverywhere a doc, doc!"
end
