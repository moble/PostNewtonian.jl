"""
    pn_order(pnsystem::PNSystem)

Return the PN order of the given `pnsystem`.

This is a `Rational{Int}` that indicates the order to which the PN expansions should be
carried out when using the given object.
"""
@export pn_order(::Type{<:PNSystem{NT,PNOrder}}) where {NT,PNOrder} = PNOrder
pn_order(pnsystem) = pn_order(typeof(pnsystem))

"""
    order_index(pnsystem::PNSystem)

Return the order index of the given `pnsystem`.

This is defined as the (one-based) index into an iterable of PN terms starting at 0pN, then
0.5pN, etc.  Specifically, this is defined as `1 + Int(2pn_order(pnsystem))`.
"""
@export order_index(pnsystem::PNSystem) = 1 + Int(2pn_order(pnsystem))

"""
    max_pn_order

The maximum PN order that can be used without overflowing the `Int` type.
"""
@export const max_pn_order = (typemax(Int) - 2) // 2

"""
    prepare_pn_order(PNOrder)

Convert the input to a half-integer of type `Rational{Int}`.

If `PNOrder` is larger than `max_pn_order`, it is set to `max_pn_order`, to avoid overflow
when computing the order index.
"""
@public function prepare_pn_order(PNOrder)
    if PNOrder < max_pn_order
        round(Int, 2PNOrder) // 2
    else
        max_pn_order
    end
end
