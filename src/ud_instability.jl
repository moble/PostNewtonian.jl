"""
    up_down_instability(pn)

Compute the range of frequencies over which the system is unstable to
increasing precession.

The returned value is a pair of dimensionless frequencies giving the lower and
upper frequencies between which we can expect instability.  If there is no
instability expected, the returned pair is just (1,1).

For compact binaries in which the spins are either aligned or anti-aligned with
the orbital angular velocity, we do not expect any precession effects — simply
by symmetry.  However, if the spin of the higher-mass object is aligned with
the orbital angular velocity and the spin of the lower-mass object is
anti-aligned, the binary is unstable to precession — meaning that any minuscule
misalignment can grow rapidly into significant precession.  This was first
reported by [Gerosa et
al. (2015)](https://link.aps.org/doi/10.1103/PhysRevLett.115.141102), and the
range over which the system is unstable is given by Eq. (2) of that reference.
We use the lowest-order approximation to convert binary separation to
frequency.  The result is also "clamped" between 0 and 1, because sometimes the
PN approximations involved break down and return unphysical values.

"""
function up_down_instability(pn)
    let ℓ̂=ℓ̂(pn.R)
        if pn.M₁ > pn.M₂
            q = pn.M₂ / pn.M₁
            χ₁ = pn.χ⃗₁ ⋅ ℓ̂
            χ₂ = pn.χ⃗₂ ⋅ ℓ̂
        else
            q = pn.M₁ / pn.M₂
            χ₂ = pn.χ⃗₁ ⋅ ℓ̂
            χ₁ = pn.χ⃗₂ ⋅ ℓ̂
        end
        if χ₁ > 0 && χ₂ < 0
            M = pn.M₁ + pn.M₂
            r₊ = M * (√(χ₁) + √(q*χ₂))^4 / (1-q)^2
            r₋ = M * (√(χ₁) - √(q*χ₂))^4 / (1-q)^2
            Ω₊ = √(M/r₊)^3
            Ω₋ = √(M/r₋)^3
            clamp.((Ω₊, Ω₋), 0, 1)
        else
            T = typeof(χ₁)
            (T(1), T(1))
        end
    end
end
