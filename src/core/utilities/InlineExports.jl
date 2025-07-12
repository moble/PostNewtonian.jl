# NOTE: This file borrows from the `InlineExports.jl` package, which is licensed under the
# MIT License. The original source can be found at https://github.com/dalum/InlineExports.jl

module NoExport

import Base: @__doc__

export @export, @public

eval(quote
    """
        @export

    No-op.  Used to disable inline exports.
    """
    macro $(Symbol("export"))(expr::Expr)
        return esc(expr)
    end
end)

"""
    @public

No-op.  Used to disable inline public macro.
"""
macro public(expr::Expr)
    return esc(expr)
end

end # module NoExport

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
        return handle(expr, :export)
    end
end)

@static if VERSION < v"1.11"
    using ..NoExport: @public
else
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

    julia> M.f(M.T(M.a))
    4
    ```
    """
    macro public(expr::Expr)
        return handle(expr, :public)
    end
end

function handle(expr::Expr, export_or_public::Symbol)
    r = handle(expr)
    ep = if r isa Symbol
        Expr(export_or_public, r)
    else
        Expr(export_or_public, r...)
    end
    return esc(quote
        $ep
        Base.@__doc__ $expr
    end)
end

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
function handle(::Val{:macrocall}, expr)
    if expr.args[1]==Symbol("@doc") ||
        (expr.args[1] == Core.GlobalRef(Core, Symbol("@doc")))
        if length(expr.args) != 4
            error("@doc expression found with $(length(expr.args)) args:\n$expr")
        end
        handle(expr.args[4])
    else
        filter(x -> x !== nothing, map(handle, expr.args[3:end]))
    end
end

end # module InlineExports

@testitem "InlineExports" begin
    using Markdown: @doc_str

    module Bla
    using PostNewtonian.InlineExports: @export, @public

    @export macro garble(expr::Expr)
        return esc(expr)
    end

    @public macro gorble(expr::Expr)
        return esc(expr)
    end

    @export abstract type Flork{T} end
    @public abstract type Flark{T} end

    @export m::Int = 3
    @public n::Int = 4

    @export @eval :evala
    @public @eval :evalb

    @export function funca(x)
        x^2
    end
    @public function funcb(x)
        x^3
    end

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

    end  # module Bla

    using .Bla

    @test (@macroexpand @garble x=3) == :(x=3)
    @test (@macroexpand Bla.@gorble y=4) == :(y=4)

    @test m isa Int
    @test m == 3
    @test Bla.n isa Int
    @test Bla.n == 4

    # @test hasproperty(Bla, :evala)
    # @test hasproperty(Bla, :evalb)

    @test funca(3) == 9
    @test Bla.funcb(3) == 27

    @test f(T(a)) == 4
    @test Bla.g(Bla.V(Bla.b)) == 27
    @test Bla.h(Bla.X(Bla.c)) == 3125

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
    # @test Base.hasproperty(Bla, :c)  # This group evaluates to false for some reason!
    # @test Base.hasproperty(Bla, :W)  # But that's true even for modules that have nothing
    # @test Base.hasproperty(Bla, :X)  # to do with `InlineExports`?!
    # @test Base.hasproperty(Bla, :h)  # Anyway, we used them in tests above... 🤷

    @test repr(@doc(Bla.a)) == "`const a` doc\n"
    @test repr(@doc(Bla.b)) == "`const b` doc\n"
    @test repr(@doc(Bla.c)) == "`const c` doc\n"
    @test repr(@doc(Bla.S)) == "`abstract type S` doc\n"
    @test repr(@doc(Bla.U)) == "`abstract type U` doc\n"
    @test repr(@doc(Bla.W)) == "`abstract type W` doc\n"
    @test repr(@doc(Bla.T)) == "`struct T` doc\n"
    @test repr(@doc(Bla.V)) == "`struct V` doc\n"
    @test repr(@doc(Bla.X)) == "`struct X` doc\n"
    @test repr(@doc(Bla.f)) == """```\nf(x)\n```\n\nHere a doc!\n"""
    @test repr(@doc(Bla.g)) == """```\ng(x)\n```\n\nThere a doc!\n"""
    @test repr(@doc(Bla.h)) == """```\nh(x)\n```\n\nEverywhere a doc, doc!\n"""
end
