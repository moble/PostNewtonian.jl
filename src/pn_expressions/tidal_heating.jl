"""
    tidal_heating(pnsystem)

Compute the rate of energy and angular-momentum absorption into each black hole in a binary.

The returned quantity is a tuple `(Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂)`, representing the corresponding
rates of change of spin (magnitude) and mass on black hole 1 and black hole 2.  These apply
to general — possibly precessing — non-eccentric binaries.  This collection of terms comes
from [Alvi (2001)](http://link.aps.org/doi/10.1103/PhysRevD.64.104020).  It probably
wouldn't be too hard to extend Alvi's analysis to eccentric systems.

Note that the validity of the result depends not only on the PN parameter ``v``, but also on
the angles of the spins relative to the separation vector ``n̂``: the smaller the angle, the
lower the ``v`` at which the approximations should be expected to break down.

See also

  - [Tagoshi and Sasaki (1994)](http://arxiv.org/abs/gr-qc/9405062)
  - [Poisson and Sasaki (1995)](https://arxiv.org/abs/gr-qc/9412027)
  - [Tagoshi et al. (1997)](https://arxiv.org/abs/gr-qc/9711072)
  - [Chatziioannou et al. (2012)](https://arxiv.org/abs/1211.1686)
"""
@compute_pn_variables function tidal_heating(pnsystem)
    # References to pages and equation numbers are from Alvi (2001)

    # Page 2, line 4
    rₕ₁ = M₁ * (1 + √(1-min(χ₁²,1)))
    rₕ₂ = M₂ * (1 + √(1-min(χ₂²,1)))

    # Page 2, line 5
    Ωₕ₁ = χ₁ / 2rₕ₁
    Ωₕ₂ = χ₂ / 2rₕ₂

    # Page 2, line 9
    b = M / v^2

    # Equivalent to Eq. (18)
    # θ is the angle between the BH's spin and the vector towards the opposite BH;
    # ϕ̇ is the rotation rate of the plane containing those two vectors.

    # Implementation note: Despite the claims of [Kahan
    # (2014)](https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf) and [Boldo
    # (2013)](https://hal.inria.fr/hal-00790071), I actually find terrible numerical
    # accuracy if I compute sin²θ using Kahan's method.  Maybe Julia is optimizing
    # something that violates the order of operations Kahan prescribed.  But I find very
    # good results using the simpler formula θ = atan √ |a⃗×b⃗|² / (a⃗⋅b⃗)², which is
    # simpler in any case, so that's what I use.

    sin²θ₁ = let cross2 = abs2vec(n̂×χ⃗₁)
        denominator = cross2 + (n̂⋅χ⃗₁)^2
        ifelse(
            iszero(denominator),
            one(cross2),
            cross2 / denominator
        )
    end
    sin²θ₂ = let cross2 = abs2vec(n̂×χ⃗₂)
        denominator = cross2 + (n̂⋅χ⃗₂)^2
        ifelse(
            iszero(denominator),
            one(cross2),
            cross2 / denominator
        )
    end
    ϕ̇₁ = let denominator = χ₁ * sin²θ₁
        ifelse(
            iszero(denominator),
            Ω,
            Ω * ℓ̂ ⋅ χ⃗₁ / denominator
        )
    end
    ϕ̇₂ = let denominator = χ₂ * sin²θ₂
        ifelse(
            iszero(denominator),
            Ω,
            Ω * ℓ̂ ⋅ χ⃗₂ / denominator
        )
    end

    # Eq. (10)
    I₀₁ = (16rₕ₁ / 5b^6) * M₁^5 * M₂^2 * sin²θ₁ * (1 - 3//4 * χ₁² + 15//4 * χ₁² * sin²θ₁)
    I₀₂ = (16rₕ₂ / 5b^6) * M₂^5 * M₁^2 * sin²θ₂ * (1 - 3//4 * χ₂² + 15//4 * χ₂² * sin²θ₂)

    # Eq. (21)
    Ṡ₁ = (ϕ̇₁ - Ωₕ₁) * I₀₁
    Ṁ₁ = ϕ̇₁ * Ṡ₁
    Ṡ₂ = (ϕ̇₂ - Ωₕ₂) * I₀₂
    Ṁ₂ = ϕ̇₂ * Ṡ₂

    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂)
end
