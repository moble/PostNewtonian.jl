@doc raw"""
Expressions from [Trestini2025](@cite)

"""
@pn_reference module Trestini2025

import PostNewtonian: G, c, v, M₁ as m₁, M₂ as m₂, M as m, r, χ₁ₗ as χ₁, χ₂ₗ as χ₂

"""
    ℵ₈(pnsystem)

4PN tidal dissipation coefficient

As given in the text below Eqs. (7.3) of [Trestini2025](@cite).  Trestini defines ℵ₂ₙ as the

> tidal dissipation (or black hole absorption) coefficient entering at the ``n``PN order
> beyond the leading, quadrupolar flux ('aleph' stands for 'absorption').
"""
@pn_expression ℵ₈(pnsystem) = 1 - 4ν + 2ν^2

# Eq. (7.3a)
@pn_expression ℱₑᴴ(pnsystem) = 32c^5 * ν^2 / 5G * ℵ₈ * y^9

end  # module Trestini2025
