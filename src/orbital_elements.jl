"""
    ℓ̂(R)

The unit vector pointing along the direction of orbital angular velocity.

"""
ℓ̂(R) = QuatVec(2R.w*R.y+2R.x*R.z, -2R.w*R.x+2R.y*R.z, R.w^2+R.z^2-R.x^2-R.y^2)
#ℓ̂(R) = QuatVec(R * 𝐤 * conj(R))


"""
    n̂(R)

The unit vector pointing from object 2 to object 1.

"""
n̂(R) = QuatVec(R.w^2+R.x^2-R.y^2-R.z^2, 2R.x*R.y+2R.w*R.z, -2R.w*R.y+2R.x*R.z)
#n̂(R) = QuatVec(R * 𝐢 * conj(R))


"""
    λ̂(R)

The unit vector pointing in the direction of the instantaneous velocity of
object 1.  This also completes the right-handed triple of (ℓ̂, n̂, λ̂).

"""
λ̂(R) = QuatVec(-2R.w*R.z+2R.x*R.y, R.w^2+R.y^2-R.x^2-R.z^2, 2R.w*R.x+2R.y*R.z)
#λ̂(R) = QuatVec(R * 𝐣 * conj(R))


Ω(;v, M=1) = v^3 / M
v(;Ω, M=1) = (M*Ω)^(1//3)
