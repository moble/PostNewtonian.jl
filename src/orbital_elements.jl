"""
    ℓ̂(R)
    ell_hat(R)

The unit vector pointing along the direction of orbital angular velocity, when
the frame is given by the rotor `R`.  This is equal to
```math
ℓ̂(R) = R ẑ R̄
```

"""
ℓ̂(R) = QuatVec(2R.w*R.y+2R.x*R.z, -2R.w*R.x+2R.y*R.z, R.w^2+R.z^2-R.x^2-R.y^2) / abs2(R)
#ℓ̂(R) = QuatVec(R * 𝐤 * conj(R))
const ell_hat = ℓ̂


"""
    n̂(R)
    n_hat(R)

The unit vector pointing from object 2 to object 1, when the frame is given by
the rotor `R`.  This is equal to
```math
n̂(R) = R x̂ R̄
```

"""
n̂(R) = QuatVec(R.w^2+R.x^2-R.y^2-R.z^2, 2R.x*R.y+2R.w*R.z, -2R.w*R.y+2R.x*R.z) / abs2(R)
#n̂(R) = QuatVec(R * 𝐢 * conj(R))
const n_hat = n̂


"""
    λ̂(R)
    lambda_hat(R)

The unit vector pointing in the direction of the instantaneous velocity of
object 1, when the frame is given by the rotor `R`.  This is equal to
```math
λ̂(R) = R ŷ R̄
```
This also completes the right-handed triple of ``(ℓ̂, n̂, λ̂)``.

"""
λ̂(R) = QuatVec(-2R.w*R.z+2R.x*R.y, R.w^2+R.y^2-R.x^2-R.z^2, 2R.w*R.x+2R.y*R.z) / abs2(R)
#λ̂(R) = QuatVec(R * 𝐣 * conj(R))
const lambda_hat = λ̂


"""
    Ω(;v, M=1)
    Omega(;v, M=1)

Orbital angular frequency.

The parameter `v` is the PN velocity parameter, and must be passed as a keyword argument — as in
`Ω(v=0.1)`.  The parameter `M` is the total mass of the binary.

"""
Ω(;v, M=1) = v^3 / M
const Omega = Ω


"""
    v(;Ω, M=1)

Post-Newtonian velocity parameter.

The parameter `Ω` is the orbital angular frequency, and must be passed as a keyword argument — as in
`v(Ω=0.1)`.  The parameter `M` is the total mass of the binary.

"""
v(;Ω, M=1) = (M*Ω)^(1//3)
