"""
    soft_fmul(a::UInt64, b::UInt64) -> UInt64

IEEE 754 double-precision multiplication on raw bit patterns.
Uses only integer operations. Bit-exact with hardware `*`.

Fully branchless: all paths computed unconditionally, results selected
via `ifelse`. Required for correct reversible circuit compilation.

The dominant cost is the 53×53 mantissa multiply, implemented as a
widening multiply via four 27×26-bit partial products (schoolbook
decomposition into half-words that fit in UInt64 without overflow).
"""
function soft_fmul(a::UInt64, b::UInt64)::UInt64
    # IEEE 754 double-precision constants
    FRAC_MASK = UInt64(0x000FFFFFFFFFFFFF)   # 52-bit stored fraction
    IMPLICIT  = UInt64(0x0010000000000000)   # bit 52 (implicit leading 1)
    SIGN_MASK = UInt64(0x8000000000000000)   # bit 63
    INF_BITS  = UInt64(0x7FF0000000000000)   # +Inf
    QNAN      = UInt64(0x7FF8000000000000)   # canonical quiet NaN
    BIAS      = Int64(1023)

    # ── Unpack ──
    sa = a >> 63
    ea = (a >> 52) & UInt64(0x7FF)
    fa = a & FRAC_MASK

    sb = b >> 63
    eb = (b >> 52) & UInt64(0x7FF)
    fb = b & FRAC_MASK

    # ── Result sign: XOR ──
    result_sign = sa ⊻ sb

    # ── Special-case predicates ──
    a_nan = (ea == UInt64(0x7FF)) & (fa != UInt64(0))
    b_nan = (eb == UInt64(0x7FF)) & (fb != UInt64(0))
    a_inf = (ea == UInt64(0x7FF)) & (fa == UInt64(0))
    b_inf = (eb == UInt64(0x7FF)) & (fb == UInt64(0))
    a_zero = (ea == UInt64(0)) & (fa == UInt64(0))
    b_zero = (eb == UInt64(0)) & (fb == UInt64(0))

    # ── Special-case results ──
    # Inf * 0 = NaN; Inf * finite = Inf; Inf * Inf = Inf
    inf_result = (result_sign << 63) | INF_BITS
    zero_result = result_sign << 63              # signed zero

    # ── Implicit leading 1 for normal, raw fraction for subnormal ──
    ma = ifelse(ea != UInt64(0), fa | IMPLICIT, fa)
    mb = ifelse(eb != UInt64(0), fb | IMPLICIT, fb)

    # Effective exponents (subnormal: stored 0 → effective 1)
    ea_eff = ifelse(ea != UInt64(0), Int64(ea), Int64(1))
    eb_eff = ifelse(eb != UInt64(0), Int64(eb), Int64(1))

    # ── Unbiased exponent sum ──
    # result_exp = ea + eb - bias
    # But we work in biased form: result_exp_biased = ea + eb - bias
    result_exp = ea_eff + eb_eff - BIAS

    # ── 53×53 → 106-bit mantissa multiply ──
    # Decompose into half-words to avoid overflow.
    # ma, mb are at most 53 bits (bit 52 = implicit 1).
    # Split each into high 27 bits and low 26 bits:
    #   ma = a_hi * 2^26 + a_lo
    #   mb = b_hi * 2^26 + b_lo
    # Product = a_hi*b_hi*2^52 + (a_hi*b_lo + a_lo*b_hi)*2^26 + a_lo*b_lo
    # Each partial product fits in 53 bits, cross terms in 54 bits.
    a_lo = ma & UInt64(0x03FFFFFF)            # low 26 bits
    a_hi = ma >> 26                           # high 27 bits
    b_lo = mb & UInt64(0x03FFFFFF)            # low 26 bits
    b_hi = mb >> 26                           # high 27 bits

    pp_ll = a_lo * b_lo                       # max 52 bits
    pp_lh = a_lo * b_hi                       # max 53 bits
    pp_hl = a_hi * b_lo                       # max 53 bits
    pp_hh = a_hi * b_hi                       # max 54 bits

    # Assemble 106-bit product as (prod_hi : prod_lo), each UInt64.
    # prod_lo holds bits [0:63], prod_hi holds bits [64:105].
    #
    # Start with pp_ll in prod_lo.
    # Add cross terms shifted left by 26.
    # Add pp_hh shifted left by 52.

    # Cross term sum (may be up to 54 bits)
    cross = pp_lh + pp_hl                     # max 54 bits, no overflow

    # prod_lo = pp_ll + (cross << 26)
    cross_lo = cross << 26                    # bits [26:79] → low 64 bits
    prod_lo = pp_ll + cross_lo
    # Carry from prod_lo addition
    carry_lo = ifelse(prod_lo < pp_ll, UInt64(1), UInt64(0))

    # prod_hi = (cross >> 38) + pp_hh + carry_lo + (pp_hh_lo_part)
    # cross >> 38 gives the high bits of (cross << 26) that didn't fit in prod_lo
    cross_hi = cross >> 38
    pp_hh_shifted = pp_hh << 52               # low bits of pp_hh*2^52 in prod_lo
    # Actually, we need to add pp_hh*2^52 to the 128-bit product.
    # pp_hh * 2^52: low 12 bits go into prod_lo[52:63], rest into prod_hi.
    # But prod_lo already has pp_ll + cross_lo. Let's redo properly.

    # Restart assembly more carefully:
    # 106-bit product P = pp_ll + (cross << 26) + (pp_hh << 52)
    #
    # Split into prod_hi (bits 64-105) and prod_lo (bits 0-63):
    #
    # Layer 1: pp_ll (52 bits) → all in prod_lo
    # Layer 2: cross << 26 → bits [26, 79]
    #   prod_lo gets bits [26,63] = cross[0:37] << 26
    #   prod_hi gets bits [64,79] = cross[38:53]
    # Layer 3: pp_hh << 52 → bits [52, 105]
    #   prod_lo gets bits [52,63] = pp_hh[0:11] << 52
    #   prod_hi gets bits [64,105] = pp_hh[12:53]

    # Let's just do it with add-with-carry.
    # acc_lo = pp_ll
    # acc_lo += cross << 26; carry1 if overflow
    # acc_lo += pp_hh << 52; carry2 if overflow
    # acc_hi = (cross >> 38) + (pp_hh >> 12) + carry1 + carry2

    acc_lo = pp_ll
    term2 = cross << 26
    sum1 = acc_lo + term2
    c1 = ifelse(sum1 < acc_lo, UInt64(1), UInt64(0))

    term3 = pp_hh << 52
    sum2 = sum1 + term3
    c2 = ifelse(sum2 < sum1, UInt64(1), UInt64(0))

    prod_lo_final = sum2
    prod_hi_final = (cross >> 38) + (pp_hh >> 12) + c1 + c2

    # ── Extract the relevant 55 bits for working format ──
    # The 106-bit product has the leading 1 at bit 105 (if both inputs had
    # implicit 1 at bit 52: 52+52+1=105, but we count from 0).
    # Actually: ma has leading 1 at bit 52, mb at bit 52.
    # Product leading 1 is at bit 104 or 105 (105 if carry).
    #
    # We need 53 mantissa bits + 3 GRS bits = 56 bits from the top of the product.
    # The product is in bits [0,105]. The top 56 bits are [50,105] or [49,105].
    #
    # Working format: bit 55 = leading 1, bits [54:3] = fraction, bits [2:0] = GRS.
    # So we want to extract starting from the leading 1, take 56 bits.
    #
    # If the product's MSB is at bit 105 (in prod_hi bit 41):
    #   wr = product[105:50] (56 bits), sticky = OR of product[49:0]
    # If MSB is at bit 104:
    #   wr = product[104:49] (56 bits), sticky = OR of product[48:0]

    # Check if bit 105 (= prod_hi bit 41) is set
    msb_at_105 = (prod_hi_final >> 41) & UInt64(1)

    # Case 1: MSB at 105 → extract bits [105:50], sticky from [49:0]
    # That's (prod_hi << 23) | (prod_lo >> 41) for the 56-bit working value... no.
    # bits [105:50] = 56 bits. prod_hi has bits [64:105] = 42 bits.
    # So bits [105:64] = prod_hi[0:41], bits [63:50] = prod_lo[63:50]
    # wr_105 = (prod_hi << 14) | (prod_lo >> 50)... let me think more carefully.
    #
    # product bit N:
    #   if N >= 64: prod_hi[N-64]
    #   else: prod_lo[N]
    #
    # We want bits [105:50] right-justified in a UInt64.
    # = (prod_hi << (64 - 50)) ... no.
    # Bits 105 down to 50. That's prod_hi[41:0] (42 bits) concatenated with
    # prod_lo[63:50] (14 bits) = 56 bits total.
    # = (prod_hi_final << 14) | (prod_lo_final >> 50)
    # But prod_hi_final might have bits above 41 — mask first.
    wr_105 = ((prod_hi_final & ((UInt64(1) << 42) - UInt64(1))) << 14) | (prod_lo_final >> 50)
    sticky_105 = ifelse((prod_lo_final & ((UInt64(1) << 50) - UInt64(1))) != UInt64(0),
                        UInt64(1), UInt64(0))
    wr_105 = wr_105 | sticky_105

    # Case 2: MSB at 104 → extract bits [104:49], sticky from [48:0]
    wr_104 = ((prod_hi_final & ((UInt64(1) << 42) - UInt64(1))) << 15) | (prod_lo_final >> 49)
    sticky_104 = ifelse((prod_lo_final & ((UInt64(1) << 49) - UInt64(1))) != UInt64(0),
                        UInt64(1), UInt64(0))
    wr_104 = wr_104 | sticky_104

    # Select and adjust exponent
    wr = ifelse(msb_at_105 != UInt64(0), wr_105, wr_104)
    # MSB at 105 means an extra factor of 2 → exponent +1
    result_exp = ifelse(msb_at_105 != UInt64(0), result_exp + Int64(1), result_exp)

    # Mask to 56 bits (bit 55 = leading 1, bits 54:0 = fraction + GRS)
    wr = wr & ((UInt64(1) << 56) - UInt64(1))

    # ── Normalize subnormal inputs ──
    # If either input was subnormal, the product's leading 1 may be lower.
    # Use the same CLZ normalization as soft_fadd.
    # The leading 1 should be at bit 55 of wr.
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

    # ── Handle subnormal result (exponent underflow) ──
    subnormal = result_exp <= Int64(0)
    shift_sub = Int64(1) - result_exp
    flush_to_zero = shift_sub >= Int64(56)
    shift_clamped = clamp(shift_sub, Int64(0), Int64(63))
    shift_u = UInt64(ifelse(flush_to_zero, Int64(0), shift_clamped))
    lost_mask_sub = (UInt64(1) << shift_u) - UInt64(1)
    lost_sub = ifelse((wr & lost_mask_sub) != UInt64(0), UInt64(1), UInt64(0))
    wr_sub_result = (wr >> shift_u) | lost_sub
    flushed_result = result_sign << 63

    wr = ifelse(subnormal,
         ifelse(flush_to_zero, wr, wr_sub_result),
         wr)
    result_exp = ifelse(subnormal, Int64(0), result_exp)

    # ── Handle overflow to ±Inf ──
    exp_overflow = result_exp >= Int64(0x7FF)
    overflow_result = (result_sign << 63) | INF_BITS

    # ── Round to nearest even ──
    guard      = (wr >> 2) & UInt64(1)
    round_bit  = (wr >> 1) & UInt64(1)
    sticky_bit = wr & UInt64(1)
    frac       = (wr >> 3) & FRAC_MASK

    grs = (guard << 2) | (round_bit << 1) | sticky_bit
    round_up = (grs > UInt64(4)) | ((grs == UInt64(4)) & ((frac & UInt64(1)) != UInt64(0)))

    frac_rounded = frac + UInt64(1)
    mant_overflow = frac_rounded == IMPLICIT
    frac_final = ifelse(round_up,
                 ifelse(mant_overflow, UInt64(0), frac_rounded),
                 frac)
    exp_after_round = ifelse(round_up & mant_overflow,
                             result_exp + Int64(1),
                             result_exp)
    exp_overflow_after_round = exp_after_round >= Int64(0x7FF)

    # ── Pack normal result ──
    exp_pack = UInt64(clamp(exp_after_round, Int64(0), Int64(0x7FE)))
    normal_result = (result_sign << 63) | (exp_pack << 52) | frac_final

    # ── Final select chain ──
    result = normal_result
    result = ifelse(exp_overflow | exp_overflow_after_round, overflow_result, result)
    result = ifelse(subnormal & flush_to_zero, flushed_result, result)
    result = ifelse(a_zero | b_zero, zero_result, result)
    result = ifelse((a_inf | b_inf) & (a_zero | b_zero), QNAN, result)  # Inf * 0
    result = ifelse((a_inf | b_inf) & !(a_zero | b_zero), inf_result, result)
    result = ifelse(a_nan | b_nan, QNAN, result)

    return result
end
