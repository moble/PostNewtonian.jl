@doc raw"""
    up_down_instability(pnsystem)

Compute the range of frequencies over which the system is unstable to increasing precession.

The returned value is a pair of dimensionless frequencies giving the lower and upper
frequencies between which we can expect instability.  If there is no instability expected,
the returned pair is just ``(1/M, 1/M)`` — where ``M`` is the total mass of the system, and
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

Note that Gerosa et al. use the convention that ``q = M_2/M_1`` — which is opposite to the
convention used in [this package](@ref q); which we account for internally in this function.
They also assume that ``M_1 \geq M_2``, which we deal with by automatically swapping the
relevant quantities.  Neither of these requires any adjustment by users of this function.
"""
@pn_expression function up_down_instability(pnsystem)
    Ωₘₐₓ = PostNewtonian.Ω(; v=c, M=M)

    # Note that, as mentioned in the docstring, Gerosa et al. use q<1
    if M₂ > M₁  # Just swap the relevant variables
        χ₂ₗ, χ₁ₗ = χ₁ₗ, χ₂ₗ
    else
        q = inv(q)
    end

    if χ₁ₗ > 0 && χ₂ₗ < 0
        # Note that `r₊ ≥ r₋`, but we keep corresponding subscripts,
        # which means that `Ω₊ ≤ Ω₋`!
        r₊ = M * (√abs(χ₁ₗ) + √abs(q * χ₂ₗ))^4 / (1 - q)^2
        r₋ = M * (√abs(χ₁ₗ) - √abs(q * χ₂ₗ))^4 / (1 - q)^2
        Ω₊ = √(M / r₊)^3
        Ω₋ = √(M / r₋)^3
        (clamp(Ω₊, zero(Ωₘₐₓ), Ωₘₐₓ), clamp(Ω₋, zero(Ωₘₐₓ), Ωₘₐₓ))
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
        v₊, v₋ = PostNewtonian.v(; Ω=Ω₊, M=M), PostNewtonian.v(; Ω=Ω₋, M=M)
        if v₁ < min(v₋, vₗᵢₘᵢₜ) && min(vₑ, vₗᵢₘᵢₜ) > v₊
            @warn (
                "\n" *
                "This system is likely to encounter the up-down instability in the frequency\n" *
                "range (Ω₊, Ω₋)=$((Ω₊, Ω₋)), which corresponds to\n" *
                "PN velocity parameters (v₊, v₋)=$((v₊, v₋)).\n" *
                "This is a true physical instability, not just a numerical issue.  Despite the\n" *
                "fact that the initial conditions contain very small precession, the system will\n" *
                "likely evolve to have very large precession.  See `up_down_instability` docs\n" *
                "for details." *
                "\n\nParameters:"
            ) M₁ M₂ χ⃗₁ χ⃗₂ R v v₁ vₑ
        end
    end
end
