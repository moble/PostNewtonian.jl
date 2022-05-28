"""
    noneccentric_evolution(m₁, m₂, χ⃗₁, χ⃗₂, Ωᵢ, Rᵢ=Rotor(true), Ω₁=Ωᵢ)

Integrate the orbital dynamics of a non-eccentric compact binary.


## Initial frequency vs. first frequency

Note the distinction between `Ωᵢ` (with subscript `i`) and `Ω₁` (with subscript
`1`).  The first, `Ωᵢ`, represents the angular velocity of the *initial
condition* from which the ODE integrator will begin; the second, `Ω₁`,
represents the target angular velocity of the first element of the output data.
That is, the ODE integration will run forwards in time from `Ωᵢ` to the merger,
and then come back to `Ωᵢ` and run backwards in time to `Ω₁`.  The output data
will stitch these two together to be one continuous (forwards-in-time) data
series.

For example, if you are trying to match to a numerical relativity (NR)
simulation, you can read the masses and spins off of the NR data when the
system is orbiting at angular velocity `Ωᵢ`.  However, you may want to know what the

Be aware that the up-down instability (where the more massive black hole has
spin aligned with the orbital angular velocity, and the less massive has spin
anti-aligned) can cause systems with nearly zero precession at the initial time
to evolve into a highly precessing system either at earlier or later times.
This is a real physical result, rather than a numerical issue.  If you want to
simulate a truly non-precessing system, you should explicitly set the in-place
components of spin to precisely 0.

"""
function noneccentric_evolution(
    M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ,
    Rᵢ=Rotor(true), Ω₁=Ωᵢ;
    PNSystem=TaylorT1
)
    vᵢ = v_Ω(Ωᵢ, M₁+M₂)
    v₁ = v_Ω(Ω₁, M₁+M₂)
    uᵢ = [  # Initial conditions for the ODE integration
        M₁;
        M₂;
        χ⃗₁.components[2:4];
        χ⃗₂.components[2:4];
        Rᵢ.components;
        vᵢ
    ]
    T = eltype(uᵢ)
    dtmin = 10√eps(T)
    reltol = √eps(T)
    abstol = √eps(T)
    Mₜₒₜ = M₁ + M₂
    Ω⃗ᵢ = QuatVec{T}(Ωᵢ * Rᵢ * imz * conj(Rᵢ))
    pn = PNSystem(uᵢ[1], uᵢ[2], QuatVec(uᵢ[3:5]...), QuatVec(uᵢ[6:8]...), Ω⃗ᵢ, vᵢ)

    estimated_time_to_merger = 5/(256ν(Mₜₒₜ) * T(vᵢ)^8) # Lowest-order PN time-to-merger
    tspan = (0estimated_time_to_merger, 4estimated_time_to_merger)
    problem_forwards = ODEProblem(noneccentric_RHS!, uᵢ, tspan, pn)

    termination_criterion = ContinuousCallback(
        (u,t,integrator) -> u[end] - 1,  # Terminate at v = 1
        terminate!,
        terminate!
    )

    solution_forwards = solve(
        problem_forwards, Tsit5(),
        reltol=reltol, abstol=abstol,
        callback=termination_criterion,
        dtmin=dtmin
    )

    if v₁ < vᵢ
        estimated_backwards_time = 5/(256ν(Mₜₒₜ) * T(v₁)^8) - estimated_time_to_merger
        tspan = (0estimated_backwards_time, -3estimated_backwards_time)
        # dtmin *= -1  # TODO: Figure out if this should happen

        problem_backwards = remake(problem_forwards; tspan=tspan)

        termination_criterion = ContinuousCallback(
            (u,t,integrator) -> u[end] - v₁,  # Terminate at v = v₁
            terminate!
        )

        solution_backwards = solve(
            problem_backwards, Tsit5(dtmin=dtmin),
            reltol=reltol, abstol=abstol,
            callback=termination_criterion
        )

        # Combine forwards and backwards
    end

    solution_forwards
end


"""
    noneccentric_RHS!(u̇, u, p, t)

Compute the right-hand side for the orbital evolution of a non-eccentric binary

Here, `u` is the ODE state vector, which can be unpacked with
[`PNDynamicalVariables`](@ref).  The parameter `p` is currently unused, but
could be used to pass un-evolved parameters through.

"""
function noneccentric_RHS!(u̇, u, pn, t)
    recalculate!(pn, u)
    u̇[1] = pn.ṁ₁
    u̇[2] = pn.ṁ₂
    u̇[3:5] = (pn.Ω⃗ᵪ₁ × QuatVec(u[3:5]...)).components[2:4]
    u̇[6:8] = (pn.Ω⃗ᵪ₂ × QuatVec(u[6:8]...)).components[2:4]
    u̇[9:12] = (pn.Ω⃗ * Quaternion(u[9:12]...) / 2).components
    u̇[13] = pn.v̇
    nothing
end
