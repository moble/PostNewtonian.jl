"""
    ℓ̂(R)

The unit vector pointing along the direction of orbital angular velocity.

"""
ℓ̂(R) = R * 𝐤 * conj(R)


"""
    n̂(R)

The unit vector pointing from object 2 to object 1.

"""
n̂(R) = R * 𝐢 * conj(R)


"""
    λ̂(R)

The unit vector pointing in the direction of the instantaneous velocity of
object 1.  This also completes the right-handed triple of (ℓ̂, n̂, λ̂).

"""
λ̂(R) = R * 𝐣 * conj(R)


Ω_v(v) = v^3
Ω_v(v, M) = Ω_v(v) / M

v_Ω(Ω::T) where T = Ω^inv(T(3))
v_Ω(Ω, M) = v_Ω(M*Ω)
