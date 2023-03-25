"""
    estimated_time_to_merger(pnsystem)

Compute the lowest-order PN approximation for the time to merger.
"""
function estimated_time_to_merger(pnsystem)
    5/(256ν(pnsystem) * v(pnsystem)^8)
end

include("dynamics/up_down_instability.jl")
include("dynamics/right_hand_sides.jl")
include("dynamics/inspiral.jl")
