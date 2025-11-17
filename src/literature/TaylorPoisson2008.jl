@doc raw"""
Expressions from [TaylorPoisson2008](@cite)

Note that there are expressions in the paper that could probably be used pretty easily to
generalize to precessing and/or eccentric systems.
"""
@pn_reference module TaylorPoisson2008

import PostNewtonian: G, c, v, M₁ as m₁, M₂ as m₂, M as m, r, χ₁ₗ as χ₁, χ₂ₗ as χ₂

# These are defined in the text below Eq. (9.5), though I've added subscripts to identify
# them as applying to the first and second bodies, respectively.
ϵ₁(pnsystem) = Base.sign(χ₁(pnsystem))
ϵ₂(pnsystem) = Base.sign(χ₂(pnsystem))

@pn_expression function Eq_9_4(pnsystem) end

@pn_expression function Eq_9_7(pnsystem) end

end  # module TaylorPoisson2008
