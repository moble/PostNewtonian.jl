let
    using Logging: Logging
    Logging.with_logger(Logging.SimpleLogger(Logging.Error)) do

        # Adapted from "Quick start" example
        M₁ = 0.4
        M₂ = 0.6
        χ⃗₁ = [0.0, 0.01, 0.1]
        χ⃗₂ = [0.01, 0.0, 0.1]
        Ωᵢ = 0.15
        Ω₁ = 99Ωᵢ / 100
        Ωₑ = 0.20
        inspiral = orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Ω₁, Ωₑ)
        uniform_in_phase(inspiral, 128)
        t′ = (inspiral.t[end] - 1.0):0.1:inspiral.t[end]
        inspiral = inspiral(t′)
        h = inertial_waveform(inspiral)

        # Exercise a couple explicit PN orders
        for PNOrder in [4//1, typemax(Int)]
            coorbital_waveform(
                orbital_evolution(M₁, M₂, χ⃗₁, χ⃗₂, Ωᵢ; Ω₁, Ωₑ, PNOrder); PNOrder
            )
        end

        # People may use the GWFrames interface; it *should* be covered above,
        # but just to be sure...
        Approximant = "TaylorT1"
        delta = (M₁ - M₂) / (M₁ + M₂)
        chi1_i = χ⃗₁
        chi2_i = χ⃗₂
        Omega_orb_i = Ωᵢ
        Omega_orb_0 = Ω₁
        GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i; Omega_orb_0)
    end
end
