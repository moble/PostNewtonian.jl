module FundamentalVariables

using ..PostNewtonian: PNState

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ

M₁(s::PNState) = @inbounds s.u[1]
M₂(s::PNState) = @inbounds s.u[2]
χ⃗₁(s::PNState) = @inbounds QuatVec(s.u[3:5]...)
χ⃗₂(s::PNState) = @inbounds QuatVec(s.u[6:8]...)
R(s::PNState) = @inbounds Rotor(s.u[9:12]...)
v(s::PNState) = @inbounds s.u[13]
Φ(s::PNState) = s.u[14]  # NO @inbounds

λ₁(::PNState) = 0
λ₂(::PNState) = 0

end
