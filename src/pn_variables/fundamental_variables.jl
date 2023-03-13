module FundamentalVariables

using ..PostNewtonian: PNSystem

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ

M₁(s::PNSystem) = @inbounds s.u[1]
M₂(s::PNSystem) = @inbounds s.u[2]
χ⃗₁(s::PNSystem) = @inbounds QuatVec(s.u[3:5]...)
χ⃗₂(s::PNSystem) = @inbounds QuatVec(s.u[6:8]...)
R(s::PNSystem) = @inbounds Rotor(s.u[9:12]...)
v(s::PNSystem) = @inbounds s.u[13]
Φ(s::PNSystem) = s.u[14]  # NO @inbounds

λ₁(::PNSystem) = 0
λ₂(::PNSystem) = 0

end
