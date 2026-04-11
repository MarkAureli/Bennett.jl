"""
    soft_fptosi(a::UInt64)::UInt64

Convert IEEE 754 double-precision float (as UInt64 bit pattern) to signed Int64
(as UInt64 bit pattern). Branchless implementation for reversible circuit compilation.

Algorithm:
1. Extract sign, exponent, mantissa from IEEE 754 encoding
2. Compute shift amount: how many bits to shift the 53-bit mantissa
3. If exponent >= 1023+52: shift left (large values)
4. If exponent < 1023+52: shift right (truncate fractional part)
5. Apply sign (two's complement negation if negative)
6. Special cases: 0 (exp=0), overflow (exp≥1086), NaN/Inf → undefined (match hardware)
"""
@inline function soft_fptosi(a::UInt64)::UInt64
    sign = (a >> 63) & UInt64(1)
    exp  = (a >> 52) & UInt64(0x7ff)
    mant = a & UInt64(0x000fffffffffffff)

    # Add implicit 1-bit for normal numbers (exp != 0)
    is_normal = ifelse(exp != UInt64(0), UInt64(1), UInt64(0))
    full_mant = mant | (is_normal << 52)  # 53-bit significand: 1.mant

    # Unbiased exponent: exp - 1023
    # Shift amount: unbiased_exp - 52 (positive = shift left, negative = shift right)
    # For the integer part, we want: full_mant >> (52 - unbiased_exp)
    #   = full_mant >> (52 - (exp - 1023))
    #   = full_mant >> (1075 - exp)
    #
    # If exp >= 1075: shift LEFT by (exp - 1075)
    # If exp < 1075: shift RIGHT by (1075 - exp)
    # If exp < 1023: result is 0 (value < 1.0)

    right_shift = UInt64(1075) - exp
    left_shift = exp - UInt64(1075)

    # Clamp shifts to valid range [0, 63]
    right_shift_clamped = ifelse(right_shift > UInt64(63), UInt64(63), right_shift)
    left_shift_clamped = ifelse(left_shift > UInt64(63), UInt64(63), left_shift)

    # Compute both paths
    result_right = full_mant >> right_shift_clamped
    result_left = full_mant << left_shift_clamped

    # Select: if exp >= 1075, use left shift; otherwise right shift
    go_left = ifelse(exp >= UInt64(1075), UInt64(1), UInt64(0))
    magnitude = ifelse(go_left == UInt64(1), result_left, result_right)

    # Zero for subnormals and zero (exp == 0 gives right_shift = 1075, so result = 0 anyway)
    # This is handled naturally by the shift.

    # Apply sign: if sign == 1, negate (two's complement: -x = ~x + 1)
    negated = (~magnitude) + UInt64(1)
    result = ifelse(sign == UInt64(1), negated, magnitude)

    return result
end
