module MathConstants

export log2, ln2, log3, ln3, log5, ln5,
    log3halves, log³╱₂, ln³╱₂, log5halves, log⁵╱₂, ln⁵╱₂,
    apery, ζ3, eulergamma, γₑ, 𝒾

Base.@irrational ln2 0.6931471805599453 log2
Base.@irrational ln3 1.0986122886681097 big"1.098612288668109691395245236922525704647490557822749451734694333637494293218608966873615754813732088787970029"
Base.@irrational ln5 1.6094379124341003 big"1.60943791243410037460075933322618763952560135426851772191264789147417898770765776463013387809317961"
Base.@irrational ln³╱₂ 0.4054651081081644 big"0.40546510810816438197801311546434913657199042346249419761401432414410067124891425126775242781731340"
Base.@irrational ln⁵╱₂ 0.9162907318741551 big"0.91629073187415506518352721176801107145010121990826246779196788198078536573796304902427055109676092"
Base.@irrational ζ3 1.2020569031595942 big"1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915"
import Base.MathConstants: eulergamma
const γₑ = eulergamma
const 𝒾 = im

# When defining @irrationals, if the third argument is a symbol, it is expected to name an
# existing constant compiled into MPFR, with prefix `mpfr_const_`, which is a very limited
# set.  All MPFR constants except `log2` already exist in `Base.MathConstants`.  Otherwise,
# that third argument should be a literal BigFloat; it cannot be a variable name.

@doc raw"""
    γₑ
    eulergamma

[Euler's constant](https://en.wikipedia.org/wiki/Euler%27s_constant) (also known as the
Euler–Mascheroni constant) is defined as the limit as ``n \to \infty`` of the difference
between the ``n``th partial sum of the harmonic series and ``\log(n)``.  This is OEIS
sequence [A001620](https://oeis.org/A001620).

```julia-repl
julia> PostNewtonian.γₑ
γ = 0.5772156649015...

julia> PostNewtonian.eulergamma
γ = 0.5772156649015...

julia> n=10_000_000; sum(1 ./ (1:n))-log(n)
0.5772157149015307
```
"""
γₑ

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
ζ3, const apery=ζ3

"""
    ln2
    log2

The natural logarithm of 2.  This is OEIS sequence
[A002162](https://oeis.org/A002162).

```julia-repl
julia> PostNewtonian.ln2
ln2 = 0.6931471805599...

julia> exp(PostNewtonian.ln2)
2.0

julia> exp(big(PostNewtonian.ln2))
2.0
```
"""
ln2, const log2=ln2

"""
    ln3
    log3

The natural logarithm of 3.  This is OEIS sequence
[A002391](https://oeis.org/A002391).

```julia-repl
julia> PostNewtonian.ln3
ln3 = 1.0986122886681...

julia> exp(PostNewtonian.ln3)
3.0000000000000004

julia> exp(big(PostNewtonian.ln3))
3.0
```
"""
ln3, const log3=ln3

"""
    ln5
    log5

The natural logarithm of 5.  This is OEIS sequence
[A016628](https://oeis.org/A016628).

```julia-repl
julia> PostNewtonian.ln5
ln5 = 1.6094379124341...

julia> exp(PostNewtonian.ln5)
4.999999999999999

julia> exp(big(PostNewtonian.ln5))
5.0
```
"""
ln5, const log5=ln5

"""
    ln³╱₂
    log³╱₂
    log3halves

The natural logarithm of 3//2.  This is OEIS sequence
[A016578](https://oeis.org/A016578).

```julia-repl
julia> PostNewtonian.ln³╱₂
ln³╱₂ = 0.4054651081081...

julia> exp(PostNewtonian.ln³╱₂)
1.5

julia> exp(big(PostNewtonian.ln³╱₂))
1.5
```
"""
ln³╱₂, const log³╱₂=ln³╱₂, const log3halves=ln³╱₂

"""
    ln⁵╱₂
    log⁵╱₂
    log5halves

The natural logarithm of 5//2.  This is OEIS sequence
[A016579](https://oeis.org/A016579).

```julia-repl
julia> PostNewtonian.ln⁵╱₂
ln⁵╱₂ = 0.9162907318741...

julia> exp(PostNewtonian.ln⁵╱₂)
2.5

julia> exp(big(PostNewtonian.ln⁵╱₂))
2.5
```
"""
ln⁵╱₂, const log⁵╱₂=ln⁵╱₂, const log5halves=ln⁵╱₂

end
