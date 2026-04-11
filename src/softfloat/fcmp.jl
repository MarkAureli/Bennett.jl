"""
    soft_fcmp_olt(a::UInt64, b::UInt64) -> UInt64

IEEE 754 ordered less-than comparison on raw bit patterns.
Returns 1 if a < b (and neither is NaN), 0 otherwise.
Fully branchless.
"""
function soft_fcmp_olt(a::UInt64, b::UInt64)::UInt64
    SIGN_MASK = UInt64(0x8000000000000000)
    ABS_MASK  = UInt64(0x7FFFFFFFFFFFFFFF)

    sa = a >> 63
    sb = b >> 63
    ea = (a >> 52) & UInt64(0x7FF)
    eb = (b >> 52) & UInt64(0x7FF)
    fa = a & FRAC_MASK
    fb = b & FRAC_MASK
    abs_a = a & ABS_MASK
    abs_b = b & ABS_MASK

    # NaN check: exponent all-ones with non-zero fraction
    a_nan = (ea == UInt64(0x7FF)) & (fa != UInt64(0))
    b_nan = (eb == UInt64(0x7FF)) & (fb != UInt64(0))
    either_nan = a_nan | b_nan

    # Both zero (±0 == ±0)
    both_zero = (abs_a == UInt64(0)) & (abs_b == UInt64(0))

    # Same sign comparison
    # For positive: a < b iff abs_a < abs_b
    # For negative: a < b iff abs_a > abs_b
    pos_lt = abs_a < abs_b      # |a| < |b|
    neg_lt = abs_a > abs_b      # |a| > |b| (more negative)

    # Different sign: negative < positive (unless both zero)
    diff_sign_lt = (sa > sb)    # a is negative, b is positive

    same_sign = sa == sb
    result = ifelse(same_sign,
                    ifelse(sa == UInt64(0), pos_lt, neg_lt),
                    diff_sign_lt)

    # Override: NaN → false, both zero → false
    result = result & (!both_zero) & (!either_nan)

    return UInt64(result)
end

"""
    soft_fcmp_oeq(a::UInt64, b::UInt64) -> UInt64

IEEE 754 ordered equal comparison on raw bit patterns.
Returns 1 if a == b (and neither is NaN), 0 otherwise.
Note: +0.0 == -0.0 per IEEE 754.
Fully branchless.
"""
function soft_fcmp_oeq(a::UInt64, b::UInt64)::UInt64
    ABS_MASK  = UInt64(0x7FFFFFFFFFFFFFFF)

    ea = (a >> 52) & UInt64(0x7FF)
    eb = (b >> 52) & UInt64(0x7FF)
    fa = a & FRAC_MASK
    fb = b & FRAC_MASK
    abs_a = a & ABS_MASK
    abs_b = b & ABS_MASK

    # NaN check
    a_nan = (ea == UInt64(0x7FF)) & (fa != UInt64(0))
    b_nan = (eb == UInt64(0x7FF)) & (fb != UInt64(0))
    either_nan = a_nan | b_nan

    # +0.0 == -0.0: both absolute values zero
    both_zero = (abs_a == UInt64(0)) & (abs_b == UInt64(0))

    # Bitwise equal or both zero
    result = (a == b) | both_zero

    # NaN != anything
    result = result & (!either_nan)

    return UInt64(result)
end

"""
    soft_fcmp_ole(a::UInt64, b::UInt64) -> UInt64

IEEE 754 ordered less-than-or-equal: a <= b and neither is NaN.
"""
@inline function soft_fcmp_ole(a::UInt64, b::UInt64)::UInt64
    return soft_fcmp_olt(a, b) | soft_fcmp_oeq(a, b)
end

"""
    soft_fcmp_une(a::UInt64, b::UInt64) -> UInt64

IEEE 754 unordered not-equal: a != b or either is NaN.
"""
@inline function soft_fcmp_une(a::UInt64, b::UInt64)::UInt64
    return UInt64(1) - soft_fcmp_oeq(a, b)
end
