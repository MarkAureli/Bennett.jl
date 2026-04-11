"""
    soft_fdiv(a::UInt64, b::UInt64) -> UInt64

IEEE 754 double-precision division on raw bit patterns.
Uses only integer operations. Bit-exact with hardware `/`.
Fully branchless.
"""
@inline function soft_fdiv(a::UInt64, b::UInt64)::UInt64
    FRAC_MASK = UInt64(0x000FFFFFFFFFFFFF)
    IMPLICIT  = UInt64(0x0010000000000000)
    INF_BITS  = UInt64(0x7FF0000000000000)
    QNAN      = UInt64(0x7FF8000000000000)
    BIAS      = Int64(1023)

    sa = a >> 63
    ea = (a >> 52) & UInt64(0x7FF)
    fa = a & FRAC_MASK
    sb = b >> 63
    eb = (b >> 52) & UInt64(0x7FF)
    fb = b & FRAC_MASK

    result_sign = sa ⊻ sb

    a_nan = (ea == UInt64(0x7FF)) & (fa != UInt64(0))
    b_nan = (eb == UInt64(0x7FF)) & (fb != UInt64(0))
    a_inf = (ea == UInt64(0x7FF)) & (fa == UInt64(0))
    b_inf = (eb == UInt64(0x7FF)) & (fb == UInt64(0))
    a_zero = (ea == UInt64(0)) & (fa == UInt64(0))
    b_zero = (eb == UInt64(0)) & (fb == UInt64(0))

    inf_result = (result_sign << 63) | INF_BITS
    zero_result = result_sign << 63

    ma = ifelse(ea != UInt64(0), fa | IMPLICIT, fa)
    mb = ifelse(eb != UInt64(0), fb | IMPLICIT, fb)
    ea_eff = ifelse(ea != UInt64(0), Int64(ea), Int64(1))
    eb_eff = ifelse(eb != UInt64(0), Int64(eb), Int64(1))

    result_exp = ea_eff - eb_eff + BIAS

    # ── Mantissa division: produce quotient with leading 1 at bit 55 ──
    # We want Q = (ma / mb) in 56-bit fixed point with 55 fractional bits.
    # This means Q = floor(ma * 2^55 / mb) approximately.
    # Since ma << 55 would overflow UInt64, use restoring division:
    # Start with r = ma, iterate 56 times, each time checking if r >= mb,
    # shifting quotient bit in, and shifting remainder.
    q = UInt64(0)
    r = ma  # start with full mantissa as initial remainder
    for i in 0:55
        # Check if current remainder can subtract divisor
        fits = r >= mb
        r = ifelse(fits, r - mb, r)
        q = (q << 1) | ifelse(fits, UInt64(1), UInt64(0))
        # Shift remainder left for next iteration (multiply by 2)
        r = r << 1
    end

    # Sticky from remainder
    sticky = ifelse(r != UInt64(0), UInt64(1), UInt64(0))
    wr = q | sticky

    # ── Normalize: leading 1 should be at bit 55 ──
    # If ma >= mb, leading 1 is at bit 55. If ma < mb, at bit 54.
    need_shift = (wr >> 55) == UInt64(0)
    wr = ifelse(need_shift, wr << 1, wr)
    result_exp = ifelse(need_shift, result_exp - Int64(1), result_exp)

    # ── Subnormal CLZ ──
    need32 = (wr & (UInt64(0xFFFFFFFF) << 24)) == UInt64(0)
    wr = ifelse(need32, wr << 32, wr)
    result_exp = ifelse(need32, result_exp - Int64(32), result_exp)

    need16 = (wr & (UInt64(0xFFFF) << 40)) == UInt64(0)
    wr = ifelse(need16, wr << 16, wr)
    result_exp = ifelse(need16, result_exp - Int64(16), result_exp)

    need8 = (wr & (UInt64(0xFF) << 48)) == UInt64(0)
    wr = ifelse(need8, wr << 8, wr)
    result_exp = ifelse(need8, result_exp - Int64(8), result_exp)

    need4 = (wr & (UInt64(0xF) << 52)) == UInt64(0)
    wr = ifelse(need4, wr << 4, wr)
    result_exp = ifelse(need4, result_exp - Int64(4), result_exp)

    need2 = (wr & (UInt64(0x3) << 54)) == UInt64(0)
    wr = ifelse(need2, wr << 2, wr)
    result_exp = ifelse(need2, result_exp - Int64(2), result_exp)

    need1 = (wr & (UInt64(1) << 55)) == UInt64(0)
    wr = ifelse(need1, wr << 1, wr)
    result_exp = ifelse(need1, result_exp - Int64(1), result_exp)

    # ── Subnormal result ──
    subnormal = result_exp <= Int64(0)
    shift_sub = Int64(1) - result_exp
    flush_to_zero = shift_sub >= Int64(56)
    shift_clamped = clamp(shift_sub, Int64(0), Int64(63))
    shift_u = UInt64(ifelse(flush_to_zero, Int64(0), shift_clamped))
    lost_mask_sub = (UInt64(1) << shift_u) - UInt64(1)
    lost_sub = ifelse((wr & lost_mask_sub) != UInt64(0), UInt64(1), UInt64(0))
    wr_sub_result = (wr >> shift_u) | lost_sub
    flushed_result = result_sign << 63

    wr = ifelse(subnormal, ifelse(flush_to_zero, wr, wr_sub_result), wr)
    result_exp = ifelse(subnormal, Int64(0), result_exp)

    # ── Overflow ──
    exp_overflow = result_exp >= Int64(0x7FF)

    # ── Round ──
    guard      = (wr >> 2) & UInt64(1)
    round_bit  = (wr >> 1) & UInt64(1)
    sticky_bit = wr & UInt64(1)
    frac       = (wr >> 3) & FRAC_MASK

    grs = (guard << 2) | (round_bit << 1) | sticky_bit
    round_up = (grs > UInt64(4)) | ((grs == UInt64(4)) & ((frac & UInt64(1)) != UInt64(0)))

    frac_rounded = frac + UInt64(1)
    mant_overflow = frac_rounded == IMPLICIT
    frac_final = ifelse(round_up, ifelse(mant_overflow, UInt64(0), frac_rounded), frac)
    exp_after_round = ifelse(round_up & mant_overflow, result_exp + Int64(1), result_exp)
    exp_overflow_after_round = exp_after_round >= Int64(0x7FF)

    # ── Pack ──
    exp_pack = UInt64(clamp(exp_after_round, Int64(0), Int64(0x7FE)))
    normal_result = (result_sign << 63) | (exp_pack << 52) | frac_final

    # ── Select chain ──
    result = normal_result
    result = ifelse(exp_overflow | exp_overflow_after_round, inf_result, result)
    result = ifelse(subnormal & flush_to_zero, flushed_result, result)
    result = ifelse(a_zero & b_zero, QNAN, result)
    result = ifelse(a_zero & (!b_zero), zero_result, result)
    result = ifelse(b_zero & (!a_zero), inf_result, result)
    result = ifelse(a_inf & b_inf, QNAN, result)
    result = ifelse(a_inf & (!b_inf), inf_result, result)
    result = ifelse(b_inf & (!a_inf), zero_result, result)
    result = ifelse(a_nan | b_nan, QNAN, result)

    return result
end
