# Declared public in src/PostNewtonian.jl
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

# Declared public in src/PostNewtonian.jl
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

@testitem "core.misc" begin
    using DoubleFloats
    using ForwardDiff: Dual
    using PostNewtonian: ζ3, γₑ, value, iscall, isadd, ismul

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
    @test iscall(:(log(1.2)), (:log,))
    @test iscall(:(log(1.2)), (:log, :sin))
    @test !iscall(:(log(1.2)), (:sin,))
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
