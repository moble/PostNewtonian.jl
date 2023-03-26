"""
    up_down_instability(pnsystem)

Compute the range of frequencies over which the system is unstable to increasing precession.

The returned value is a pair of dimensionless frequencies giving the lower and upper
frequencies between which we can expect instability.  If there is no instability expected,
the returned pair is just (1/M, 1/M) — where ``M`` is the total mass of the system, and
``1/M`` is the upper limit of physically reasonable frequencies.

For compact binaries in which the spins are either aligned or anti-aligned with the orbital
angular velocity, we do not expect any precession effects — simply by symmetry.  However, if
the spin of the higher-mass object is aligned with the orbital angular velocity and the spin
of the lower-mass object is anti-aligned, the binary is unstable to precession — meaning
that any minuscule misalignment can grow rapidly into significant precession.  This was
first reported by [Gerosa et al. (2015)](http://arxiv.org/abs/1506.09116), and the range
over which the system is unstable is given by Eq. (2) of that reference.  We use the
lowest-order approximation to convert binary separation to frequency.  The result is also
"clamped" between ``0`` and ``1/M``, because sometimes the PN approximations involved break
down and return values outside of those physically plausible limits.
"""
@pn_expression function up_down_instability(pnsystem)
    T = eltype(pnsystem)
    Ωₘₐₓ = PostNewtonian.Ω(v=1, M=M)
    if M₂ ≤ M₁
        q, χ₂ₗ, χ₁ₗ = inv(q), χ₁ₗ, χ₂ₗ
    end
    if χ₁ₗ > 0 && χ₂ₗ < 0
        r₊ = M * (√(χ₁ₗ) + √abs(q*χ₂ₗ))^4 / (1-q)^2
        r₋ = M * (√(χ₁ₗ) - √abs(q*χ₂ₗ))^4 / (1-q)^2
        Ω₊ = √(M/r₊)^3
        Ω₋ = √(M/r₋)^3
        clamp.((Ω₊, Ω₋), T(0), Ωₘₐₓ)
    else
        (Ωₘₐₓ, Ωₘₐₓ)
    end
end


@doc raw"""
    function up_down_instability_warn(pnsystem, v₁, vₑ, vₗᵢₘᵢₜ=1//2)

If this system is likely to encounter the up-down instability, log a warning with details.

This function issues the warning if the system is reasonably non-precessing (``\chi_\perp
\leq 10^{-2}``) in its current configuration (as given by `pnsystem`), and the range of
frequencies `(v₁, vₑ)` over which it will be integrated is likely to encounter the up-down
instability — except that frequencies above `vₗᵢₘᵢₜ` will be ignored, as PN is likely to
have broken down anyway.

See [`up_down_instability`](@ref) for details of the calculation of the unstable region.
"""
@pn_expression function up_down_instability_warn(pnsystem, v₁, vₑ, vₗᵢₘᵢₜ=1//2)
    if χₚₑᵣₚ ≤ 1e-2 && !iszero(χₚₑᵣₚ)
        (Ω₊, Ω₋) = up_down_instability(pnsystem)
        v₊, v₋ = PostNewtonian.v(Ω=Ω₊, M=M), PostNewtonian.v(Ω=Ω₋, M=M)
        if v₁ < min(v₋, vₗᵢₘᵢₜ) && min(vₑ, vₗᵢₘᵢₜ) > v₊
            @warn (
                "This system is likely to encounter the up-down instability in the\n"
                * "frequency range (Ω₊, Ω₋)=$((Ω₊, Ω₋)),\ncorresponding "
                * "to the range of PN velocity parameters (v₊, v₋)=$((v₊, v₋)).\n"
                * "This is a true physical instability; not just a numerical issue.\n"
                * "Despite the initial conditions containing very small precession,\n"
                * "the system will likely evolve to have very large precession."
            ) M₁ M₂ χ⃗₁ χ⃗₂ R v v₁ vₑ
        end
    end
end
