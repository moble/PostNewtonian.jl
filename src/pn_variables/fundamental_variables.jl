module FundamentalVariables

using ..PostNewtonian: AbstractPNSystem

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ

M₁(s::AbstractPNSystem) = @inbounds s.u[1]
M₂(s::AbstractPNSystem) = @inbounds s.u[2]
χ⃗₁(s::AbstractPNSystem) = @inbounds QuatVec(s.u[3:5]...)
χ⃗₂(s::AbstractPNSystem) = @inbounds QuatVec(s.u[6:8]...)
R(s::AbstractPNSystem) = @inbounds Rotor(s.u[9:12]...)
v(s::AbstractPNSystem) = @inbounds s.u[13]
Φ(s::AbstractPNSystem) = s.u[14]  # NO @inbounds

λ₁(::AbstractPNSystem) = 0
λ₂(::AbstractPNSystem) = 0

end
