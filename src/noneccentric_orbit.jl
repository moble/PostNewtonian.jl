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
    PNSys=TaylorT1, PNOrder=7//2,
    check_up_down_instability=true,
)
    if Ω₁ > Ωᵢ
        error(
            "Initial frequency Ωᵢ=$Ωᵢ should be greater"
            * " than or equal to first frequency Ω₁=$Ω₁"
        )
    end

    vᵢ = v(Ω=Ωᵢ, M=M₁+M₂)
    v₁ = v(Ω=Ω₁, M=M₁+M₂)
    uᵢ = [  # Initial conditions for the ODE integration
        M₁;
        M₂;
        χ⃗₁.vec;
        χ⃗₂.vec;
        Rᵢ.components;
        vᵢ
    ]
    T = eltype(uᵢ)
    dtmin = 10√eps(T)
    reltol = √eps(T)
    abstol = √eps(T)
    Mₜₒₜ = M₁ + M₂
    Ω⃗ᵢ = QuatVec{T}(Ωᵢ * Rᵢ * imz * conj(Rᵢ))
    pn = PNSys(PNOrder, T)
    unpack!(pn, uᵢ)

    if check_up_down_instability
        χₚₑᵣₚ = let n̂=n̂(pn.R),λ̂=λ̂(pn.R)
            √((pn.χ⃗₁ ⋅ n̂)^2 + (pn.χ⃗₁ ⋅ λ̂)^2 + (pn.χ⃗₂ ⋅ n̂)^2 + (pn.χ⃗₂ ⋅ λ̂)^2)
        end
        if χₚₑᵣₚ ≤ 1e-2
            (Ω₊, Ω₋) = up_down_instability(pn)
            if Ω₁ < Ω₋ < 1//8 || Ω₁ < Ω₊ < 1//8
                @warn (
                    "This system is likely to encounter the up-down instability in the\n"
                    * "frequency range (Ω₊, Ω₋)=$((Ω₊, Ω₋)).\n"
                    * "This is a true physical instability; not just a numerical issue.\n"
                    * "Despite the initial conditions containing very small precession,\n"
                    * "the system will likely evolve to have very large precession."
                )
            end
        end
    end

    estimated_time_to_merger = 5/(256ν(M₁, M₂) * T(vᵢ)^8) # Lowest-order PN time-to-merger
    tspan = (T(0), 4estimated_time_to_merger)
    problem_forwards = ODEProblem(noneccentric_RHS!, uᵢ, tspan, pn)

    termination_criterion = ContinuousCallback(
        (u,t,integrator) -> u[end] - 1,  # Terminate at v = 1
        terminate!
    )

    solution_forwards = solve(
        problem_forwards, AutoVern9(Rodas5()),
        reltol=reltol, abstol=abstol,
        callback=termination_criterion,
        dtmin=dtmin
    )

    if v₁ < vᵢ
        estimated_backwards_time = 5/(256ν(M₁, M₂) * T(v₁)^8) - estimated_time_to_merger
        tspan = (T(0), -3estimated_backwards_time)
        # dtmin *= -1  # TODO: Figure out if this should happen

        problem_backwards = remake(problem_forwards; tspan=tspan)

        termination_criterion = ContinuousCallback(
            (u,t,integrator) -> u[end] - v₁,  # Terminate at v = v₁
            terminate!
        )

        solution_backwards = solve(
            problem_backwards, AutoVern9(Rodas5()),
            reltol=reltol, abstol=abstol,
            callback=termination_criterion,
            dtmin=dtmin
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
    u̇[1] = pn.Ṁ₁
    u̇[2] = pn.Ṁ₂
    u̇[3:5] = pn.χ⃗̇₁.vec
    u̇[6:8] = pn.χ⃗̇₂.vec
    u̇[9:12] = pn.Ṙ.components
    u̇[13] = pn.v̇
    nothing
end
