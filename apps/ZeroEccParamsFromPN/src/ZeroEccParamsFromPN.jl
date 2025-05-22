module ZeroEccParamsFromPN

using PostNewtonian: PostNewtonian
import Quaternionic: QuatVecF64
using ArgParse: ArgParse

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
  - If d₀ is given, figure out the corresponding Ω₀ and proceed as below
  - If tₘ is given, search for the value of Ω₀ that gives the right time
    and compute the corresponding d₀, ȧ₀, and Nₒ
  - If Nₒ is given, search for the value of Ω₀ that gives the right number
    of orbits and compute the corresponding d₀, ȧ₀, and tₘ
  - If Ω₀ is given, compute the corresponding d₀, ȧ₀, tₘ, and Nₒ
* If Dᵣ is given, compute Ωᵣ and proceed as below
* If Ωᵣ is given, evolve the given system forwards in time, then
  - If d₀ is given, compute the corresponding Ω₀, ȧ₀, tₘ, and Nₒ

"""
function julia_main(args=nothing)::Cint
    # @info "Args:" args
    try
        # Parse command line arguments; arguments not provided are set to `nothing`
        parsed_args = parse_commandline(args)
        q::Float64OrNothing = parsed_args["q"]
        χ⃗₁::QuatVecF64OrNothing = parsed_args["chiA"]
        χ⃗₂::QuatVecF64OrNothing = parsed_args["chiB"]
        Ω₀::Float64OrNothing = parsed_args["Omega0"]
        d₀::Float64OrNothing = parsed_args["D0"]
        tₘ::Float64OrNothing = parsed_args["tMerger"]
        Nₒ::Float64OrNothing = parsed_args["NOrbits"]
        Ωᵣ::Float64OrNothing = parsed_args["OmegaRef"]
        Dᵣ::Float64OrNothing = parsed_args["DRef"]

        # @info "Parsed:" parsed_args

        M₁ = q/(1+q)
        M₂ = 1/(1+q)

        zero_ecc_params_from_pn(M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, d₀, tₘ, Nₒ, Ωᵣ, Dᵣ)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function zero_ecc_params_from_pn(M₁, M₂, χ⃗₁, χ⃗₂, Ω₀, d₀, tₘ, Nₒ, Ωᵣ, Dᵣ)
    # End every integration at Ω=0.1; this comes from the original script,
    # spec/Support/Python/ZeroEccParamsFromPN.py, and is kept here for consistency; there is
    # no "right" choice.
    vₑ = PostNewtonian.v(; Ω=0.1) # ≈ 0.464

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

    # Handle input D0 because this easily transforms to the problem with input Omega0
    if !isnothing(d₀)
        v₀ = PostNewtonian.separation_inverse(d₀, pnsystem)
        pnsystem.state[PostNewtonian.vindex] = v₀
        Ω₀ = PostNewtonian.Ω(pnsystem)
    end

    # r′₀ is a gauge choice; try a few values and see what the difference is
    r′₀s = (1.0, 10.0)  # These are values from the old script

    # TODO: Incorporate r′₀
    r′₀ = r′₀s[1]

    if isnothing(Dᵣ) && isnothing(Ωᵣ)
        if !isnothing(Ω₀)
            r₀ = PostNewtonian.r(pnsystem)
            ȧ₀ = PostNewtonian.ṙ(pnsystem) / r₀

            pnevolution = PostNewtonian.orbital_evolution(pnsystem; vₑ)
            Nₒ = pnevolution[:Φ, end] / 2π
            tₘ = pnevolution.t[end]

            println("###############################")
            println("Results for rPrime0 = $(r′₀):")
            println("Omega0 = $Ω₀")
            println("D0 = $r₀")
            println("adot0 = $ȧ₀")
            println("Approximate nOrbits = $Nₒ")
            println("Approximate tMerger = $tₘ")
        elseif !isnothing(Nₒ)
            throw(ErrorException("Not implemented: Nₒ"))
            # # First find an omega0 that gives the right number of orbits, then use that to
            # # find d₀, adot0.
            # function helperFunc(args)
            #     omega0 = args[1]
            #     println("nOrbits not implemented")
            #     return abs(nOrbits(q, χ⃗₁, χ⃗₂, omega0) - Nₒ)
            # end
            # omega = fmin(helperFunc, 0.01)[1]
            # # fromΩ₀(omega, r′₀s, q, χ⃗₁, χ⃗₂)
            # println("fromΩ₀ not implemented")
        elseif !isnothing(tₘ)
            throw(ErrorException("Not implemented: tₘ"))
            # # First find an omega0 that gives the right time, then use that to find d₀,
            # # adot0.
            # function helperFunc(args)
            #     omega0 = args[1]
            #     println("totalTime not implemented")
            #     return abs(totalTime(q, χ⃗₁, χ⃗₂, omega0) - tₘ)
            # end
            # omega = fmin(helperFunc, 0.01)[1]
            # # fromΩ₀(omega, r′₀s, q, χ⃗₁, χ⃗₂)
            # println("fromΩ₀ not implemented")
        end
    else
        throw(ErrorException("Not implemented: Dᵣ / Ωᵣ"))
    end
end

export main
(@main)(args) = julia_main(args)

end # module ZeroEccParamsFromPN
