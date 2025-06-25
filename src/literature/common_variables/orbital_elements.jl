"""
    n̂(pnsystem)
    n̂(R)
    n_hat(pnsystem)
    n_hat(R)

The unit vector pointing from object 2 to object 1, when the frame is given by the rotor
`R`.  This is equal to
```math
n̂(R) = R x̂ R̄
```
"""
n̂(R) = QuatVec(R(𝐢))
n̂(s::VecOrPNSystem) = n̂(R(s))
const n_hat = n̂

"""
    λ̂(pnsystem)
    λ̂(R)
    lambda_hat(pnsystem)
    lambda_hat(R)

The unit vector pointing in the direction of the instantaneous velocity of object 1, when
the frame is given by the rotor `R`.  This is equal to
```math
λ̂(R) = R ŷ R̄
```
This also completes the right-handed triple of ``(n̂, λ̂, ℓ̂)``.
"""
λ̂(R) = QuatVec(R(𝐣))
λ̂(s::VecOrPNSystem) = λ̂(R(s))
const lambda_hat = λ̂

"""
    ℓ̂(pnsystem)
    ℓ̂(R)
    ell_hat(pnsystem)
    ell_hat(R)

The unit vector pointing along the direction of orbital angular velocity, when the frame is
given by the rotor `R`.  This is equal to
```math
ℓ̂(R) = R ẑ R̄
```
"""
ℓ̂(R) = QuatVec(R(𝐤))
ℓ̂(s::VecOrPNSystem) = ℓ̂(R(s))
const ell_hat = ℓ̂

@doc raw"""
    Ω(pnsystem)
    Ω(;v, M=1)
    Omega(pnsystem)
    Omega(;v, M=1)

Orbital angular frequency.

The parameter `v` is the PN velocity parameter, and must be passed as a keyword argument —
as in `Ω(v=0.1)`.  The parameter `M` is the total mass of the binary.  They are related *by
definition* as
```math
\Omega \colonequals \frac{v^3}{M}.
```
See also [`v`](@ref).
"""
Ω(; v, M=1) = v^3 / M
Ω(s::VecOrPNSystem) = Ω(; v=v(s), M=M(s))
const Omega = Ω

lnv(s::VecOrPNSystem) = ln(v(s))
