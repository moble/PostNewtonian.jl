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
        help = "Orbital angular frequency at which q, chiA, chiB have the given values."
        arg_type = Float64
        required = false
        "--DRef"
        help = "Separation distance at which q, chiA, chiB have the given values."
        arg_type = Float64
        required = false
    end

    if isnothing(args)
        return ArgParse.parse_args(s)
    else
        return ArgParse.parse_args(args, s)
    end
end

const Float64OrNothing = Union{Float64,Nothing}
const QuatVecF64OrNothing = Union{QuatVecF64,Nothing}

"""
* If Dᵣ and Ωᵣ are not given
  - If r₀ is given, figure out the corresponding Ω₀ and proceed as below
  - If tₘ is given, search for the value of Ω₀ that gives the right time and compute the
    corresponding r₀, ȧ₀, and Nₒ
  - If Nₒ is given, search for the value of Ω₀ that gives the right number of orbits and
    compute the corresponding r₀, ȧ₀, and tₘ
  - If Ω₀ is given, compute the corresponding r₀, ȧ₀, tₘ, and Nₒ
* If Dᵣ is given, compute Ωᵣ and proceed as below
* If Ωᵣ is given, evolve the given system forwards in time, then
  - If r₀ is given, compute the corresponding Ω₀, ȧ₀, tₘ, and Nₒ


Test this using something like one of the following commands, run from the
apps/ZeroEccParamsFromPN directory:
```bash
julia --project -e 'using ZeroEccParamsFromPN' -- --q=4.3 --chiA=0.1,0.2,0.3 --chiB=0.3,0.2,0.1 --Omega0=0.01
julia --project -e 'using ZeroEccParamsFromPN' -- --q=4.3 --chiA=0.1,0.2,0.3 --chiB=0.3,0.2,0.1 --D0=20
julia --project -e 'using ZeroEccParamsFromPN' -- --q=4.3 --chiA=0.1,0.2,0.3 --chiB=0.3,0.2,0.1 --tMerger=10000
julia --project -e 'using ZeroEccParamsFromPN' -- --q=4.3 --chiA=0.1,0.2,0.3 --chiB=0.3,0.2,0.1 --NOrbits=22
```
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
        M₁ = q/(1+q)
        M₂ = 1/(1+q)

        zero_ecc_params_from_pn(M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, r₀, tₘ, Nₒ, Ωᵣ, Dᵣ)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function zero_ecc_params_from_pn(M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, r₀, tₘ, Nₒ, Ωᵣ, Dᵣ)
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
            pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
            while pnevolution[:Φ, end] / 2π < Nₒ
                v₀ *= 0.9
                pnsystem.state[PostNewtonian.vindex] = v₀
                pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
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
                    pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
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
            # as the unbracketed one, so if the latter gets fixed, we could try to switch
            # algorithms.
            Ω₀ = find_zero(
                Ω₀ -> begin
                    pnsystem.state[PostNewtonian.vindex] = PostNewtonian.v(; Ω=Ω₀)
                    pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
                    tₘ - pnevolution.t[end]
                end,
                (Ω₀/2, 0.9Ωₑ),
            )
        end

        if !isnothing(r₀)
            for r′₀ ∈ r′₀s
                v₀ = PostNewtonian.r⁻¹(r₀, pnsystem, r′₀)
                pnsystem.state[PostNewtonian.vindex] = v₀
                Ω₀ = PostNewtonian.Ω(pnsystem)
                evolve_and_evaluate(Ω₀, pnsystem, r′₀, vₑ)
            end
        elseif !isnothing(Ω₀)
            for r′₀ ∈ r′₀s
                evolve_and_evaluate(Ω₀, pnsystem, r′₀, vₑ)
            end
        else
            # This is an error.  I don't see how this could happen, but just in case...
            throw(ErrorException("Ω₀ has not been given or calculated"))
        end
    else
        throw(ErrorException("Not implemented: Dᵣ / Ωᵣ"))
    end
end

function evolve_and_evaluate(Ω₀, pnsystem, r′₀, vₑ)
    pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
    Nₒ = pnevolution[:Φ, end] / 2π
    tₘ = pnevolution.t[end]
    r₀ = PostNewtonian.r(pnsystem, r′₀)  # This is redundant if r₀ is given, but that's fine
    ȧ₀ = PostNewtonian.ṙ(pnsystem, r′₀) / r₀

    println("###############################")
    println("Results for rPrime0 = $(r′₀):")
    println("Omega0 = $Ω₀")
    println("D0 = $r₀")
    println("adot0 = $ȧ₀")
    println("Approximate nOrbits = $Nₒ")
    println("Approximate tMerger = $tₘ")
end

export main
(@main)(args) = julia_main(args)

end # module ZeroEccParamsFromPN
