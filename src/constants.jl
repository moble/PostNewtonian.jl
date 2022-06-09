import Base.MathConstants: eulergamma

using Base: @irrational
# When defining @irrationals, if the third argument is a symbol, it is expected
# to name an existing constant compiled into MPFR, with prefix `mpfr_const_`,
# which is a very limited set.  All such constants except `log2` already exist
# in `Base.MathConstants`.  Otherwise, it should be a BigFloat; it cannot be a
# variable name.
Base.@irrational log2 0.6931471805599453 log2
Base.@irrational log3 1.0986122886681097 big"1.098612288668109691395245236922525704647490557822749451734694333637494293218608966873615754813732088787970029"
Base.@irrational log5 1.6094379124341003 big"1.60943791243410037460075933322618763952560135426851772191264789147417898770765776463013387809317961"
Base.@irrational log3halves 0.4054651081081644 big"0.40546510810816438197801311546434913657199042346249419761401432414410067124891425126775242781731340"
Base.@irrational log5halves 0.9162907318741551 big"0.91629073187415506518352721176801107145010121990826246779196788198078536573796304902427055109676092"
Base.@irrational apery 1.2020569031595942 big"1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915"

"""
    ζ3
    apery

[Apéry's constant](https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_constant) is
defined as ``ζ(3)``, where ``ζ`` is the Riemann zeta function.  It is OEIS
sequence [A002117](https://oeis.org/A002117).

"""
apery, const ζ3=apery
