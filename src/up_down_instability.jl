"""
    up_down_instability(u)
    up_down_instability(M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v)

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
reported by [Gerosa et al. (2015)](http://arxiv.org/abs/1506.09116), and the
range over which the system is unstable is given by Eq. (2) of that reference.
We use the lowest-order approximation to convert binary separation to
frequency.  The result is also "clamped" between 0 and 1, because sometimes the
PN approximations involved break down and return unphysical values.

"""
function up_down_instability(M₁, M₂, χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ, χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ, Rʷ, Rˣ, Rʸ, Rᶻ, v)
    χ⃗₁ = QuatVec(χ⃗₁ˣ, χ⃗₁ʸ, χ⃗₁ᶻ)
    χ⃗₂ = QuatVec(χ⃗₂ˣ, χ⃗₂ʸ, χ⃗₂ᶻ)
    R = Quaternion(Rʷ, Rˣ, Rʸ, Rᶻ)
    let ℓ̂=ℓ̂(R)
        if M₁ > M₂
            q = M₂ / M₁
            χ₁ = χ⃗₁ ⋅ ℓ̂
            χ₂ = χ⃗₂ ⋅ ℓ̂
        else
            q = M₁ / M₂
            χ₂ = χ⃗₁ ⋅ ℓ̂
            χ₁ = χ⃗₂ ⋅ ℓ̂
        end
        if χ₁ > 0 && χ₂ < 0
            M = M₁ + M₂
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
up_down_instability(u) = up_down_instability(u...)
