"""
    FDPNSystem(asdf)

TODO UPDATE

"""
function FDPNSystem(::Type{PN}, PNOrder) where {NT,ST,PN<:PNSystem{NT,ST}}
    raw_PN = parameterless_type(PN)
    return raw_PN{FastDifferentiation.Node,
                  DenseVector{FastDifferentiation.Node},
                  prepare_pn_order(PNOrder)}(
                      [FastDifferentiation.Node(s) for s ∈ symbols(raw_PN)])
end

function FDPNSystem(t::Type{PN}) where {NT,ST,PNOrder,PN<:PNSystem{NT,ST,PNOrder}}
    return FDPNSystem(t, PNOrder)
end

function FDPNSystem(pn::PN) where {PN<:PNSystem}
    return FDPNSystem(typeof(pn))
end
