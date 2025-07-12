@irrational ζ3 1.2020569031595942 big"1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915"
"""
    ζ3
    apery

[Apéry's constant](https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_constant) is
defined as ``ζ(3)``, where ``ζ`` is the Riemann zeta function.  This is OEIS
sequence [A002117](https://oeis.org/A002117).

```julia-repl
julia> PostNewtonian.apery
ζ3 = 1.2020569031595...

julia> PostNewtonian.ζ3
ζ3 = 1.2020569031595...

julia> sum((1:10_000_000).^-3)
1.2020569031595896
```
"""
ζ3  # We document it this way because `@irrational` cannot handle docstrings.

@doc raw"""
    γₑ

[Euler's constant](https://en.wikipedia.org/wiki/Euler%27s_constant) (also known as the
Euler–Mascheroni constant) is defined as the limit as ``n \to \infty`` of the difference
between the ``n``th partial sum of the harmonic series and ``\log(n)``.  This is OEIS
sequence [A001620](https://oeis.org/A001620).

This is distinct from the Euler's *number* ``e``, which is defined as the limit as ``n \to
\infty`` of the sum of ``1/n!``.

Note that it is usually denoted simply as `γ` in the broader literature (and in Julia's own
`Base.MathConstants`), but that symbol is used in the post-Newtonian literature for the
quantity denoted in this package as [`γₚₙ`](@ref).  To distinguish between the two, the PN
literature uses ``\gamma_\mathrm{E}`` for Euler's constant.  (There is no Unicode subscript
"E", so we use "e" instead.)

```julia-repl
julia> PostNewtonian.γₑ
γₑ = 0.5772156649015...

julia> n=10_000_000; sum(1 ./ (1:n))-log(n)
0.5772157149015307
```
"""
γₑ  # We document it here because docs don't fit nicely in the main file where it's defined.

"""
    value(x)

Return `x` or the value wrapped by the `Dual` number `x`
"""
@public value(x::T) where {T} = hasfield(T, :value) ? getfield(x, :value) : x

"""
    find_symbols_of_type(mod, T)

Given a module `mod` (not just its name, but the actual imported module), find all objects
inside that module that are instances of the given type `T`.  The returned quantity is a
vector of `Symbol`s.
"""
function find_symbols_of_type(mod, T)
    return filter(n -> getproperty(mod, n) isa T, names(mod))
end

"""
    iscall(x, symbols)

Return `true` if the `Expr` `x` is a call to any element of `symbols`.
"""
iscall(x, symbols) = MacroTools.isexpr(x, :call) && x.args[1] ∈ symbols

"""
    isadd(x)

Return `true` if the `Expr` `x` is a call to `(+)` or `:+`.
"""
isadd(x) = iscall(x, ((+), :+))

"""
    ismul(x)

Return `true` if the `Expr` `x` is a call to `(*)` or `:*`.
"""
ismul(x) = iscall(x, ((*), :*))

"""
    flatten_binary!(expr, symbols)

Flatten nested binary operations — that is, apply associativity repeatedly.
"""
function flatten_binary!(expr, symbols)
    while iscall(expr, symbols) && any(x -> iscall(x, symbols), expr.args[2:end])
        args = expr.args[2:end]
        i₊ = findfirst(x -> iscall(x, symbols), args)
        args′ = [first(symbols); args[1:(i₊ - 1)]; args[i₊].args[2:end]; args[(i₊ + 1):end]]
        expr.args[:] = args′[1:length(expr.args)]
        append!(expr.args, args′[(1 + length(expr.args)):end])
    end
    return expr
end

flatten_add!(expr) = flatten_binary!(expr, ((+), :+))
flatten_mul!(expr) = flatten_binary!(expr, ((*), :*))

"""
    apply_to_first_add!(expr, func)

Apply `func` to the first sub-expression found in a "prewalk"-traversal of `expr` that
satisfies [`isadd`](@ref).  If `func` acts in place, so does this function.  In either case,
the expression should be returned.
"""
function apply_to_first_add!(expr, func)
    found_add = false
    MacroTools.prewalk(expr) do x
        if !found_add && isadd(x)
            found_add = true
            func(x)
        else
            x
        end
    end
end

@testitem "core.misc" begin
    using DoubleFloats
    using ForwardDiff: Dual
    using PostNewtonian: ζ3, γₑ, value, find_symbols_of_type
    using PostNewtonian:
        iscall, isadd, ismul, flatten_add!, flatten_mul!, apply_to_first_add!

    # ζ3 formula
    @test ζ3 ≈ sum((1:10_000_000) .^ -3)

    # γₑ formula
    N = 10_000_000
    δγₑ = 1/2N
    @test γₑ ≈ sum(1 ./ (1:N)) - log(N) - δγₑ

    # value
    struct Dummy
        value
    end
    struct Dummier
        valuet
    end
    for T ∈ (Float16, Float32, Float64, Double16, Double32, Double64, BigFloat)
        x = T(big"1.2")
        @test value(x) isa T
        @test value(x) == x
        d = Dummy(x)
        @test value(d) isa T
        @test value(d) == x
        dr = Dummier(x)
        @test value(dr) isa Dummier
        @test value(dr) == dr
        v = (T(big"3.4"), T(big"5.6"), T(big"7.8"))
        for i ∈ eachindex(v)
            ẋ = Dual{:Taggo}(x, v[begin:i])
            @test value(ẋ) isa T
            @test value(ẋ) == x
        end
    end

    # iscall and friends
    @test isadd(:(1 + 2))
    @test isadd(:(1 + 2c + 3c^2))
    @test isadd(:(1 + (2n + 3m)))
    @test isadd(Expr(:call, :+, 1, 2))
    @test isadd(Expr(:call, (+), 1, 2))
    @test !isadd(:(1 * 2))
    @test ismul(:(1 * 2))
    @test ismul(:(1 * 2c * 3c^2))
    @test ismul(:(1 * (2c * 3d)))
    @test ismul(Expr(:call, :*, 1, 2))
    @test ismul(Expr(:call, (*), 1, 2))
    @test !ismul(:(1 + 2))
end
