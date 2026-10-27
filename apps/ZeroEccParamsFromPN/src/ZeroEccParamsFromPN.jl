module ZeroEccParamsFromPN

import PostNewtonian: PostNewtonian
import Quaternionic: QuatVecF64
import ArgParse: ArgParse
import Roots: find_zero
import DataInterpolations: CubicSpline

function ArgParse.parse_item(::Type{QuatVecF64}, x::AbstractString)
    components = split(x, ",")
    if length(components) != 3
        throw(ArgumentError("Expected three-component vector, got $(length(components))"))
    end
    try
        return QuatVecF64(
            parse(Float64, components[1]),
            parse(Float64, components[2]),
            parse(Float64, components[3]),
        )
    catch e
        throw(ArgumentError("Failed to parse components of x=$x as `Float64`s:\n$e"))
    end
end

function parse_commandline(args=nothing)
    s = ArgParse.ArgParseSettings(;
        description="Find low-eccentricity initial parameters (D0, Omega0, adot0) with PN.",
        version="$(pkgversion(PostNewtonian))",
    )

    ArgParse.add_arg_group!(
        s, "required arguments (pass all of the following)"; required=true
    )
    ArgParse.@add_arg_table! s begin
        "--q"
        help = "Mass ratio MA/MB"
        arg_type = Float64
        required = true
        "--chiA"
        help = "Spin vector of black hole A (comma-separated)"
        arg_type = QuatVecF64
        required = true
        "--chiB"
        help = "Spin vector of black hole B (comma-separated)"
        arg_type = QuatVecF64
        required = true
    end

    ArgParse.add_arg_group!(
        s,
        "criterion arguments (pass exactly one of the following)";
        exclusive=true,
        required=true,
    )
    ArgParse.@add_arg_table! s begin
        "--D0"
        help = "Initial separation distance"
        arg_type = Float64
        "--Omega0"
        help = "Initial orbital angular frequency"
        arg_type = Float64
        "--tMerger"
        help = "Time until merger — approximate; could take a minute."
        arg_type = Float64
        "--NOrbits"
        help = "Number of orbits until merger — approximate; could take a minute."
        arg_type = Float64
    end

    ArgParse.add_arg_group!(
        s,
        "reference criteria (pass at most one of the following)";
        exclusive=true,
        required=false,
    )
    ArgParse.@add_arg_table! s begin
        "--OmegaRef"
        help = (
            "[Not yet implemented] Orbital angular frequency at which q, chiA, chiB have " *
            "the given values."
        )
        arg_type = Float64
        required = false
        "--DRef"
        help = (
            "[Not yet implemented] Separation distance at which q, chiA, chiB have the " *
            "given values."
        )
        arg_type = Float64
        required = false
    end

    # Add a single optional argument --quiet with default value false
    ArgParse.add_arg_group!(s, "optional behavior argument"; required=false)
    ArgParse.@add_arg_table! s begin
        "--skip_check_up_down_instability"
        help = "Don't check for the up-down instability."
        action = :store_true
        "--experimental_D0_hack"
        help = (
            "Add experimental fit for PN-NR offset to D0 (subject to change): `v₀^5*χₑ + " *
            "(1/8*ν*χₑ + 10/3)*v₀^3 + 3/7*v₀*ν - 1/33`"
        )
        action = :store_true
    end

    if isnothing(args)
        return ArgParse.parse_args(s)
    else
        return ArgParse.parse_args(args, s)
    end
end

const Float64OrNothing = Union{Float64,Nothing}
const BoolOrNothing = Union{Bool,Nothing}

"""
This is the main function that will be called when the script is run.
"""
function julia_main(args=nothing)::Cint
    try
        # Parse command line arguments; arguments not provided are set to `nothing`
        parsed_args = parse_commandline(args)
        q::Float64 = parsed_args["q"]
        χ⃗₁::QuatVecF64 = parsed_args["chiA"]
        χ⃗₂::QuatVecF64 = parsed_args["chiB"]
        Ω₀::Float64OrNothing = parsed_args["Omega0"]
        r₀::Float64OrNothing = parsed_args["D0"]
        tₘ::Float64OrNothing = parsed_args["tMerger"]
        Nₒ::Float64OrNothing = parsed_args["NOrbits"]
        Ωᵣ::Float64OrNothing = parsed_args["OmegaRef"]
        Dᵣ::Float64OrNothing = parsed_args["DRef"]
        skipud::BoolOrNothing = parsed_args["skip_check_up_down_instability"]
        D0_hack::BoolOrNothing = parsed_args["experimental_D0_hack"]
        M₁ = q/(1+q)
        M₂ = 1/(1+q)

        zero_ecc_params_from_pn(M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, r₀, tₘ, Nₒ, Ωᵣ, Dᵣ, skipud, D0_hack)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function zero_ecc_params_from_pn(
    M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, r₀, tₘ, Nₒ, Ωᵣ, Dᵣ, skipud, D0_hack, quiet=false
)
    if all(isnothing, (Ω₀, r₀, tₘ, Nₒ))
        throw(ArgumentError("At least one of Ω₀, r₀, tₘ, or Nₒ must be given"))
    end

    check_up_down_instability = !skipud

    # End every integration at Ω=0.1; this comes from the original script,
    # spec/Support/Python/ZeroEccParamsFromPN.py, and is kept here for consistency; there is
    # no "right" choice.
    Ωₑ = 0.1
    vₑ = PostNewtonian.v(; Ω=Ωₑ) # ≈ 0.464

    # If Ω₀ is not given, we just choose a sensible default that will be long enough that we
    # don't typically need to backtrack, but not so long that it will be slow to integrate.
    v₀ = if isnothing(Ω₀)
        # Using estimated_time_to_merger = 5M / (256ν * v^8), we solve for v assuming a very
        # long time to merger (which would probably never be simulated in practice).
        long_time_to_merger = 200_000  # M
        ν = PostNewtonian.ν(M₁, M₂)
        (5 / (256ν * long_time_to_merger))^(1//8) # ≈ 0.16 — 0.20 for q ≲ 20
    else
        PostNewtonian.v(; Ω=Ω₀)
    end

    # Note that this system may have the wrong v₀
    pnsystem = PostNewtonian.BBH(; M₁, M₂, χ⃗₁, χ⃗₂, v=v₀)

    # r′₀ is a gauge choice; try a few values and see what the difference is
    r′₀s = (1.0, 10.0)  # These are values from the old script

    if isnothing(Dᵣ) && isnothing(Ωᵣ)
        if !isnothing(Nₒ)
            # One way to get a first guess is to just evolve with the naive value of `Ω₀`,
            # then find the `t` corresponding to the `Φ` that gives the right number of
            # orbits.  Then run `find_zero` to find the value of `Ω₀` that gives the right
            # number of orbits.

            # Evolve naive system
            pnevolution = PostNewtonian.orbital_evolution(
                pnsystem; vₑ, check_up_down_instability
            )
            while pnevolution[:Φ, end] / 2π < Nₒ
                v₀ *= 0.9
                pnsystem.state[PostNewtonian.vindex] = v₀
                pnevolution = PostNewtonian.orbital_evolution(
                    pnsystem; vₑ, check_up_down_instability
                )
            end
            Nₑ = pnevolution[:Φ, end] / 2π
            spline = CubicSpline(pnevolution[:v], pnevolution[:Φ] / 2π)
            v₀ = spline(Nₑ - Nₒ)
            pnsystem.state[PostNewtonian.vindex] = v₀
            Ω₀ = PostNewtonian.Ω(pnsystem)

            # Now actively search for the value of Ω₀ that gives the right Nₒ.
            Ω₀ = find_zero(
                Ω₀ -> begin
                    pnsystem.state[PostNewtonian.vindex] = PostNewtonian.v(; Ω=Ω₀)
                    pnevolution = PostNewtonian.orbital_evolution(
                        pnsystem; vₑ, check_up_down_instability
                    )
                    Nₑ = pnevolution[:Φ, end] / 2π
                    Nₒ - Nₑ
                end,
                (Ω₀/2, 0.9Ωₑ),
            )

        elseif !isnothing(tₘ)
            # Establish a rough first guess
            v₀ = (5 / (256PostNewtonian.ν(M₁, M₂) * tₘ))^(1//8)
            pnsystem.state[PostNewtonian.vindex] = v₀
            Ω₀ = PostNewtonian.Ω(pnsystem)
            # Now actively search for the value of Ω₀ that gives the right tₘ.
            #
            # TODO: We have to use a bracketed search because it is too easy for the
            # algorithm to choose a value of Ω₀ larger than Ωₑ, which means that the
            # terminator never triggers.  The bracketing algorithm is typically not as fast
            # as the unbracketed one, so if issue #87 gets fixed, we could try to switch
            # algorithms.
            Ω₀ = find_zero(
                Ω₀ -> begin
                    pnsystem.state[PostNewtonian.vindex] = PostNewtonian.v(; Ω=Ω₀)
                    pnevolution = PostNewtonian.orbital_evolution(
                        pnsystem; vₑ, check_up_down_instability
                    )
                    tₘ - pnevolution.t[end]
                end,
                (Ω₀/2, 0.9Ωₑ),
            )
        end

        if !isnothing(r₀)
            return [
                convert_evolve_and_evaluate(
                    r₀, pnsystem, r′₀, vₑ, check_up_down_instability, D0_hack, quiet
                ) for r′₀ ∈ r′₀s
            ]
        elseif !isnothing(Ω₀)
            return [
                evolve_and_evaluate(
                    Ω₀, pnsystem, r′₀, vₑ, check_up_down_instability, D0_hack, quiet
                ) for r′₀ ∈ r′₀s
            ]
        else
            # This is an error.  I don't see how this could happen, but just in case...
            throw(ErrorException("Ω₀ has not been given or calculated"))
        end
    else
        throw(ErrorException("Not implemented: Dᵣ / Ωᵣ"))
    end
end

function convert_evolve_and_evaluate(
    r₀, pnsystem, r′₀, vₑ, check_up_down_instability, D0_hack, quiet
)
    v₀ = if D0_hack
        # The input r₀ is really r₀ᴺᴿ, but we want to find r₀ᴾᴺ according to
        #
        #    γ₀ᴾᴺ ≈ γ₀ᴺᴿ + δγ
        #
        # and thence the value of Ω₀ᴾᴺ with which to evolve the system.  By construction,
        # v₀ᴾᴺ=v₀ᴺᴿ, but we don't know how to calculate v from γ₀ᴺᴿ; we need γ₀ᴾᴺ to
        # calculate it.  But we need v to calculate γ₀ᴾᴺ, so we have to use fixed-point
        # iteration.
        M = PostNewtonian.M(pnsystem)
        ν = PostNewtonian.ν(pnsystem)
        χₑ = PostNewtonian.χₑ(pnsystem)
        γ₀ᴺᴿ = M / r₀
        γ₀ᴾᴺ = γ₀ᴺᴿ  # Initial guess for γ₀ᴾᴺ
        for i ∈ 1:100  # Had-code limit for number of iterations to avoid infinite loops
            v = PostNewtonian.γₚₙ⁻¹(γ₀ᴾᴺ, pnsystem, r′₀)
            δγ = v^2 * (v^5 * χₑ + (1/8 * ν * χₑ + 10/3) * v^3 + 3/7 * v * ν - 1/33)
            γ₀ᴾᴺ′ = γ₀ᴺᴿ + δγ
            if abs(γ₀ᴾᴺ′ - γ₀ᴾᴺ) < 10 * eps(γ₀ᴾᴺ)
                #@info "Breaking at iteration $i with" γ₀ᴾᴺ γ₀ᴺᴿ δγ
                break  # Convergence criterion
            end
            γ₀ᴾᴺ = γ₀ᴾᴺ′
            #@info "Iteration $i" γ₀ᴾᴺ γ₀ᴺᴿ δγ
        end
        PostNewtonian.γₚₙ⁻¹(γ₀ᴾᴺ, pnsystem, r′₀)
    else
        PostNewtonian.r⁻¹(r₀, pnsystem, r′₀)
    end
    pnsystem.state[PostNewtonian.vindex] = v₀

    # let # Debug correction
    #     v = PostNewtonian.v(pnsystem)
    #     M = PostNewtonian.M(pnsystem)
    #     ν = PostNewtonian.ν(pnsystem)
    #     χₑ = PostNewtonian.χₑ(pnsystem)
    #     δγ = v^2 * (v^5 * χₑ + (1/8 * ν * χₑ + 10/3) * v^3 + 3/7 * v * ν - 1/33)
    #     γᴾᴺ = PostNewtonian.γₚₙ(pnsystem, r′₀)
    #     γᴺᴿ = γᴾᴺ - δγ
    #     @info "Debugging correction" r₀ M/γᴺᴿ PostNewtonian.r(pnsystem, r′₀)
    # end

    Ω₀ = PostNewtonian.Ω(pnsystem)
    return evolve_and_evaluate(
        Ω₀, pnsystem, r′₀, vₑ, check_up_down_instability, D0_hack, quiet
    )
end

function evolve_and_evaluate(
    Ω₀, pnsystem, r′₀, vₑ, check_up_down_instability, D0_hack, quiet
)
    pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ, check_up_down_instability)
    Nₒ = pnevolution[:Φ, end] / 2π
    tₘ = pnevolution.t[end]

    r₀, ȧ₀ = if D0_hack
        # We want the output D0 to correspond to the NR value, but we have calculated the
        # PN value.  We have a correction term such that
        #     γ₀ᴺᴿ ≈ γ₀ᴾᴺ - v₀^2 * (...)
        # where γ=M/r is the inverse of the separation distance.
        ν = PostNewtonian.ν(pnsystem)
        M = PostNewtonian.M(pnsystem)
        χₑ = PostNewtonian.χₑ(pnsystem)
        γ₀ᴾᴺ = PostNewtonian.γₚₙ(pnsystem, r′₀)
        v₀ = PostNewtonian.v(; Ω=Ω₀)
        v̇ = -PostNewtonian.𝓕(pnsystem) / PostNewtonian.𝓔′(pnsystem)
        γ₀ᴺᴿ =
            γ₀ᴾᴺ - v₀^2 * (v₀^5 * χₑ + (1/8 * ν * χₑ + 10/3) * v₀^3 + 3/7 * v₀ * ν - 1/33)
        r₀ᴺᴿ = M / γ₀ᴺᴿ
        # γ = M/r, so γ̇ = -M ṙ₀ / r₀² = -ṙ₀ γ₀² / M, and ṙ = -M γ̇ / γ₀².
        # We can write ȧ = ṙ γ / M = -γ̇/γ.
        γ̇₀ᴾᴺ = PostNewtonian.γ̇ₚₙ(pnsystem)
        γ̇₀ᴺᴿ =
            γ̇₀ᴾᴺ -
            v̇ *
            v₀ *
            (7 * v₀^5 * χₑ + (1/8 * ν * χₑ + 10/3) * 5 * v₀^3 + 9/7 * v₀ * ν - 2/33)
        ȧ₀ᴺᴿ = -γ̇₀ᴺᴿ / γ₀ᴺᴿ
        r₀ᴺᴿ, ȧ₀ᴺᴿ
    else
        r₀ = PostNewtonian.r(pnsystem, r′₀)  # This is redundant if r₀ is given; that's fine
        ȧ₀ = PostNewtonian.ṙ(pnsystem, r′₀) / r₀
        r₀, ȧ₀
    end

    if !quiet
        println("###############################")
        println("Results for rPrime0 = $(r′₀):")
        println("Omega0 = $Ω₀")
        println("D0 = $r₀")
        println("adot0 = $ȧ₀")
        println("Approximate nOrbits = $Nₒ")
        println("Approximate tMerger = $tₘ")
    end

    return (r′₀, Ω₀, r₀, ȧ₀, Nₒ, tₘ)
end

export main
(@main)(args) = julia_main(args)

end # module ZeroEccParamsFromPN
