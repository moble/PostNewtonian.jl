module FundamentalVariables

using ..PostNewtonian: PNState

export M₁, M₂, χ⃗₁, χ⃗₂, R, v

M₁(s::PNState) = s.u[1]
M₂(s::PNState) = s.u[2]
χ⃗₁(s::PNState) = QuatVec(s.u[3:5]...)
χ⃗₂(s::PNState) = QuatVec(s.u[6:8]...)
R(s::PNState) = Rotor(s.u[9:12]...)
v(s::PNState) = s.u[13]

end
