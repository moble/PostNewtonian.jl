@doc raw"""
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


# Notes

The hardest parts to compute here involve the relative angles between the spins and the
black-hole separation vectors.  Alvi constructs a spherical coordinate system centered on
each black hole where the ``z`` axis is given by the direction of the spin, and ``\theta``
and ``\phi`` represent the direction to the other black hole.  While he makes a (somewhat
ambiguous) choice about the origin of the ``\phi`` coordinate, only ``\dot{\phi}`` comes
into the equations, so we don't really care about that origin.

Note that Alvi uses ``\mathbf{n}`` to represent the "normal to the orbital plane", whereas
we — and most of the rest of the post-Newtonian literature — use [``\hat{\ell}``](@ref
PostNewtonian.ℓ̂) for this vector and [``\hat{n}``](@ref PostNewtonian.n̂) to represent the
separation vector pointing *from* object 2 *to* object 1.  For convenience, we define
```math
\hat{n}_i = \begin{cases}
-\hat{n} & \text{i=1}, \\
\hphantom{-}\hat{n} & \text{i=2}.
\end{cases}
```

Alvi's construction is also somewhat adiabatic, so we treat the spins and orbital plane as
constant in the calculation of the instantaneous tidal heating — though they evolve slowly
over time — and the angles ``\theta`` and ``\phi`` as evolving rapidly.  In our formulation,
then, the goal is to find the angle ``\theta_i(t)`` between ``\vec{\chi}_i`` and
``\hat{n}_i``, and the rotation rate ``\dot{\phi}_i(t)`` of ``\hat{n}_i`` about the
``\vec{\chi}_i`` axis.

Computing the angle between vectors is a somewhat infamously tricky problem.  There are
various [claims](https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf) floating around
about the best ways to compute quantities involving areas and angles of triangles.  While
these claims are [surely true for areas](https://inria.hal.science/hal-00790071), I am more
skeptical of the relevance for angles.  I find it best to realize that you probably don't
need the angle *per se*, but trigonometric functions of the angle — like ``\sin^2
\theta_i``, which is what we actually need in this case.  In particular, I believe the best
results come from computing
```math
\sin^2\theta_i = \frac{
    \left|\hat{n}_i \times \vec{\chi}_i\right|^2
}{
    \left| \vec{\chi}_i\right|^2
}.
```
If the denominator is zero, we set ``\sin^2\theta_i = 1`` from physical considerations.

Now consider the quantity ``\hat{n}_i \times \vec{\chi}_i``.  We next aim to calculate the
rotation rate ``\dot{\phi}_i`` of this vector about ``\vec{\chi}_i``.  We begin by directly
calculating
```math
\partial_t \left(\hat{n}_i \times \vec{\chi}_i\right)
=
\left(\partial_t \hat{n}_i\right) \times \vec{\chi}_i
=
\left(\Omega\, \hat{\ell} \times \hat{n}_i\right) \times \vec{\chi}_i
=
\mp \Omega\, \hat{\lambda} \times \vec{\chi}_i,
```
where the negative sign is chosen for ``i=1`` and positive for ``i=2``.  Now, from more
fundamental considerations, we can understand the components of this change.  Since we
assume that ``\vec{\chi}_i`` is constant at each instant for the purposes of calculation
here, the only way for ``\hat{n}_i \times \vec{\chi}_i`` to change is either because
``\hat{n}_i`` rotates about ``\vec{\chi}_i``, or because the angle ``\theta_i`` is changing.
We can express this as
```math
\partial_t \left(\hat{n}_i \times \vec{\chi}_i\right)
=
\left( \dot{\phi}\, \hat{\chi}_i \right) \times \left(\hat{n}_i \times \vec{\chi}_i\right)
+
\dot{\theta_i}\, \cot \theta_i\, \left(\hat{n}_i \times \vec{\chi}_i\right).
```
Since these two components are orthogonal, we can obtain ``\dot{\phi}`` directly by taking
the component of this quantity along ``\hat{\chi}_i \times \left(\hat{n}_i \times
\vec{\chi}_i\right)``:
```math
\dot{\phi}_i
=
\frac{
    \left( \Omega\, \hat{\lambda} \times \vec{\chi}_i \right)
    \cdot
    \left(
        \vphantom{\hat{\lambda}} \hat{\chi}_i
        \times
        \left(\hat{n} \times \vec{\chi}_i\right)
    \right)
}{
    \left| \hat{\chi}_i \times \left(\hat{n} \times \vec{\chi}_i\right) \right|^2
}
=
\Omega\, \frac{\hat{\ell} \cdot \hat{\chi}_i}{\sin^2 \theta_i}.
```
Here again, we may run into numerical trouble if ``\left| \vec{\chi}_1 \right| \approx 0``,
in which case we again use physical arguments to take ``\dot{\phi}_i = \Omega``.  We might
also expect to run into trouble if ``\sin^2 \theta_i \approx 0``, which corresponds to a
polar orbit, in which case Alvi's approximations break down.  This turns out to not be a
problem *numerically*, because of the cancellation with the numerator, except when ``\sin^2
\theta_i = 0`` exactly.  In this case, we choose ``\dot{\phi}_i = 0``.

Note that the sign of ``\hat{n}_i`` has dropped out of the calculations of both
``\sin^2\theta_i`` and ``\dot{\phi}_i``, cancelling with the signs that had appeared next to
``\Omega``.
"""
@pn_expression function tidal_heating(pnsystem)
    # References to pages and equation numbers are from Alvi (2001)

    # Eq. (21)
    Ṡ₁ = Î₀₁ * v^12 * @pn_expansion (v^3/M * ϕ̇̂₁ - Ωₕ₁)
    Ṁ₁ = Î₀₁ * v^12 * v^3/M * ϕ̇̂₁ * @pn_expansion 5 (v^3/M * ϕ̇̂₁ - Ωₕ₁)  # ϕ̇₁ * Ṡ₁
    Ṡ₂ = Î₀₂ * v^12 * @pn_expansion (v^3/M * ϕ̇̂₂ - Ωₕ₂)
    Ṁ₂ = Î₀₂ * v^12 * v^3/M * ϕ̇̂₂ * @pn_expansion 5 (v^3/M * ϕ̇̂₂ - Ωₕ₂)  # ϕ̇₂ * Ṡ₂

    (Ṡ₁, Ṁ₁, Ṡ₂, Ṁ₂)
end
