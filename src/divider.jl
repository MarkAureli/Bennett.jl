"""
    soft_udiv(a::UInt64, b::UInt64) -> UInt64

Unsigned integer division via restoring division algorithm.
Uses only operations our pipeline already handles (add, sub, shift, compare, select).
Fully branchless. Returns a ÷ b (truncated).
"""
function soft_udiv(a::UInt64, b::UInt64)::UInt64
    q = UInt64(0)
    r = UInt64(0)
    for i in 63:-1:0
        # Shift remainder left, bring in bit i of a
        r = (r << 1) | ((a >> i) & UInt64(1))
        # Trial: can we subtract b?
        fits = r >= b
        r = ifelse(fits, r - b, r)
        q = ifelse(fits, q | (UInt64(1) << i), q)
    end
    return q
end

"""
    soft_urem(a::UInt64, b::UInt64) -> UInt64

Unsigned integer remainder. Returns a % b.
"""
function soft_urem(a::UInt64, b::UInt64)::UInt64
    r = UInt64(0)
    for i in 63:-1:0
        r = (r << 1) | ((a >> i) & UInt64(1))
        fits = r >= b
        r = ifelse(fits, r - b, r)
    end
    return r
end
