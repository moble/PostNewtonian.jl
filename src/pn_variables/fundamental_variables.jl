module FundamentalVariables

using ..PostNewtonian: AbstractPNSystem

export M₁, M₂, χ⃗₁, χ⃗₂, R, v, Φ

M₁(s::AbstractPNSystem) = @inbounds s.state[1]
M₂(s::AbstractPNSystem) = @inbounds s.state[2]
χ⃗₁(s::AbstractPNSystem) = @inbounds QuatVec(s.state[3:5]...)
χ⃗₂(s::AbstractPNSystem) = @inbounds QuatVec(s.state[6:8]...)
R(s::AbstractPNSystem) = @inbounds Rotor(s.state[9:12]...)
v(s::AbstractPNSystem) = @inbounds s.state[13]
Φ(s::AbstractPNSystem) = s.state[14]  # NO @inbounds

λ₁(::AbstractPNSystem) = 0
λ₂(::AbstractPNSystem) = 0

end
