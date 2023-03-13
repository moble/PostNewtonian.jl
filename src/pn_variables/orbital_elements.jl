"""
    n̂(R)
    n_hat(R)

The unit vector pointing from object 2 to object 1, when the frame is given by the rotor
`R`.  This is equal to
```math
n̂(R) = R x̂ R̄
```
"""
n̂(R) = QuatVec(R(𝐢))
n̂(s::PNSystem) = n̂(R(s))
const n_hat = n̂


"""
    λ̂(R)
    lambda_hat(R)

The unit vector pointing in the direction of the instantaneous velocity of object 1, when
the frame is given by the rotor `R`.  This is equal to
```math
λ̂(R) = R ŷ R̄
```
This also completes the right-handed triple of ``(n̂, λ̂, ℓ̂)``.
"""
λ̂(R) = QuatVec(R(𝐣))
λ̂(s::PNSystem) = λ̂(R(s))
const lambda_hat = λ̂


"""
    ℓ̂(R)
    ell_hat(R)

The unit vector pointing along the direction of orbital angular velocity, when the frame is
given by the rotor `R`.  This is equal to
```math
ℓ̂(R) = R ẑ R̄
```
"""
ℓ̂(R) = QuatVec(R(𝐤))
ℓ̂(s::PNSystem) = ℓ̂(R(s))
const ell_hat = ℓ̂


"""
    Ω(;v, M=1)
    Omega(;v, M=1)

Orbital angular frequency.

The parameter `v` is the PN velocity parameter, and must be passed as a keyword argument —
as in `Ω(v=0.1)`.  The parameter `M` is the total mass of the binary.
"""
Ω(;v, M=1) = v^3 / M
Ω(s::PNSystem) = Ω(v=v(s), M=M(s))
const Omega = Ω


"""
    v(;Ω, M=1)

Post-Newtonian velocity parameter.

The parameter `Ω` is the orbital angular frequency, and must be passed as a keyword argument
— as in `v(Ω=0.1)`.  The parameter `M` is the total mass of the binary.
"""
v(;Ω, M=1) = (M*Ω)^(1//3)


lnv(s::PNSystem) = ln(v(s))
