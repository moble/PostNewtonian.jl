"""
    up_down_instability(pnsystem)

Compute the range of frequencies over which the system is unstable to increasing precession.

The returned value is a pair of dimensionless frequencies giving the lower and upper
frequencies between which we can expect instability.  If there is no instability expected,
the returned pair is just (1,1).

For compact binaries in which the spins are either aligned or anti-aligned with the orbital
angular velocity, we do not expect any precession effects — simply by symmetry.  However, if
the spin of the higher-mass object is aligned with the orbital angular velocity and the spin
of the lower-mass object is anti-aligned, the binary is unstable to precession — meaning
that any minuscule misalignment can grow rapidly into significant precession.  This was
first reported by [Gerosa et al. (2015)](http://arxiv.org/abs/1506.09116), and the range
over which the system is unstable is given by Eq. (2) of that reference.  We use the
lowest-order approximation to convert binary separation to frequency.  The result is also
"clamped" between 0 and 1, because sometimes the PN approximations involved break down and
return unphysical values.
"""
@compute_pn_variables function up_down_instability(pnsystem)
    T = typeof(χ₁ₗ)
    if M₂ ≤ M₁
        q = inv(q)
        χ₂ₗ, χ₁ₗ = χ₁ₗ, χ₂ₗ
    end
    if χ₁ₗ > 0 && χ₂ₗ < 0
        r₊ = M * (√(χ₁ₗ) + √abs(q*χ₂ₗ))^4 / (1-q)^2
        r₋ = M * (√(χ₁ₗ) - √abs(q*χ₂ₗ))^4 / (1-q)^2
        Ω₊ = √(M/r₊)^3
        Ω₋ = √(M/r₋)^3
        clamp.((Ω₊, Ω₋), T(0), T(1))
    else
        (T(1), T(1))
    end
end
