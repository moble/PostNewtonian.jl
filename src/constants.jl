import Base.MathConstants: eulergamma

# When defining @irrationals, if the third argument is a symbol, it is expected
# to name an existing constant compiled into MPFR, with prefix `mpfr_const_`,
# which is a very limited set.  All such constants except `log2` already exist
# in `Base.MathConstants`.
Base.@irrational log2 0.6931471805599453 log2
Base.@irrational apery 1.2020569031595942 big"1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915"

"""
    ζ3
    apery

[Apéry's constant](https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_constant) is
defined as ``ζ(3)``, where ``ζ`` is the Riemann zeta function.  It is OEIS
sequence [A002117](https://oeis.org/A002117).

"""
apery, const ζ3=apery
