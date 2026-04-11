"""
    soft_fadd(a::UInt64, b::UInt64) -> UInt64

IEEE 754 double-precision addition on raw bit patterns.
Uses only integer operations. Bit-exact with hardware `+`.

Fully branchless: every path is computed unconditionally and results are
selected via `ifelse`. This ensures LLVM emits `select` instructions
(not `br`+`phi`), which is required for correct reversible circuit
compilation — branching causes false-path sensitization bugs in the
phi resolution algorithm.

Working format: mantissa shifted left by 3 for guard/round/sticky bits.
  bit 55 = implicit 1, bits 54-3 = 52-bit fraction, bits 2/1/0 = G/R/S.
"""
function soft_fadd(a::UInt64, b::UInt64)::UInt64
    # IEEE 754 double-precision constants
    FRAC_MASK = UInt64(0x000FFFFFFFFFFFFF)   # 52-bit stored fraction
    IMPLICIT  = UInt64(0x0010000000000000)   # bit 52 (implicit leading 1)
    SIGN_MASK = UInt64(0x8000000000000000)   # bit 63
    INF_BITS  = UInt64(0x7FF0000000000000)   # +Inf
    QNAN      = UInt64(0x7FF8000000000000)   # canonical quiet NaN

    # ── Unpack ──
    sa = a >> 63
    ea = (a >> 52) & UInt64(0x7FF)
    fa = a & FRAC_MASK

    sb = b >> 63
    eb = (b >> 52) & UInt64(0x7FF)
    fb = b & FRAC_MASK

    # ── Special-case predicates (computed unconditionally) ──
    a_nan = (ea == UInt64(0x7FF)) & (fa != UInt64(0))
    b_nan = (eb == UInt64(0x7FF)) & (fb != UInt64(0))
    a_inf = (ea == UInt64(0x7FF)) & (fa == UInt64(0))
    b_inf = (eb == UInt64(0x7FF)) & (fb == UInt64(0))
    a_zero = (ea == UInt64(0)) & (fa == UInt64(0))
    b_zero = (eb == UInt64(0)) & (fb == UInt64(0))

    # ── Special-case results ──
    # NaN: any NaN input → QNAN
    # Inf+Inf same sign → Inf, Inf-Inf → NaN
    inf_inf_result = ifelse(sa == sb, a, QNAN)
    # Inf + finite → Inf
    inf_finite_a = a   # a is Inf
    inf_finite_b = b   # b is Inf
    # Zero + Zero: same sign → keep sign, diff sign → +0
    zero_zero_result = ifelse(sa == sb, a, UInt64(0))

    # ── Order by magnitude: ensure |a| >= |b| ──
    a_mag = a & ~SIGN_MASK
    b_mag = b & ~SIGN_MASK
    swap = a_mag < b_mag

    sa_ord = ifelse(swap, sb, sa)
    sb_ord = ifelse(swap, sa, sb)
    ea_ord = ifelse(swap, eb, ea)
    eb_ord = ifelse(swap, ea, eb)
    fa_ord = ifelse(swap, fb, fa)
    fb_ord = ifelse(swap, fa, fb)

    # ── Implicit leading 1 for normal numbers ──
    ma = ifelse(ea_ord != UInt64(0), fa_ord | IMPLICIT, fa_ord)
    mb = ifelse(eb_ord != UInt64(0), fb_ord | IMPLICIT, fb_ord)

    # Effective exponents (subnormal: stored 0 → effective 1)
    ea_eff = ifelse(ea_ord != UInt64(0), ea_ord, UInt64(1))
    eb_eff = ifelse(eb_ord != UInt64(0), eb_ord, UInt64(1))
    d = ea_eff - eb_eff                        # >= 0 since |a| >= |b|

    # ── Working format: shift left by 3 for G/R/S room ──
    wa = ma << 3
    wb = mb << 3

    # ── Align: shift wb right by d, tracking sticky ──
    # Case 1: d >= 56 — everything shifts out, just sticky
    wb_large = ifelse(wb != UInt64(0), UInt64(1), UInt64(0))
    # Case 2: 0 < d < 56 — shift with sticky tracking
    # Clamp d to valid shift range for computing the mask (avoid UB at d=0 or d>=64)
    d_clamped = ifelse(d == UInt64(0), UInt64(1), ifelse(d >= UInt64(64), UInt64(63), d))
    lost_mask = (UInt64(1) << d_clamped) - UInt64(1)
    sticky = ifelse((wb & lost_mask) != UInt64(0), UInt64(1), UInt64(0))
    wb_mid = (wb >> d) | sticky
    # Case 3: d == 0 — no shift
    wb_aligned = ifelse(d >= UInt64(56), wb_large,
                 ifelse(d > UInt64(0),   wb_mid,
                                         wb))

    # ── Add AND subtract mantissas (both computed unconditionally) ──
    wr_add = wa + wb_aligned
    wr_sub = wa - wb_aligned

    same_sign = sa_ord == sb_ord
    wr_raw = ifelse(same_sign, wr_add, wr_sub)

    # Exact cancellation: subtraction yields zero → result is +0.0
    exact_cancel = (!same_sign) & (wr_sub == UInt64(0))

    result_sign = sa_ord

    # Use wr_raw for normalization (if exact_cancel, we'll override at the end)
    # Substitute 1 for zero to avoid undefined normalization behavior
    wr = ifelse(exact_cancel, UInt64(1), wr_raw)
    result_exp = Int64(ea_eff)

    # ── Normalize: overflow (bit 56 set from addition carry) ──
    overflow = (wr >> 56) != UInt64(0)
    lost_ov = wr & UInt64(1)
    wr_ov = (wr >> 1) | lost_ov                # preserve sticky in bit 0
    exp_ov = result_exp + Int64(1)
    wr = ifelse(overflow, wr_ov, wr)
    result_exp = ifelse(overflow, exp_ov, result_exp)

    # ── Normalize: underflow (leading 1 below bit 55) ──
    # Six constant-shift stages (binary-search CLZ)
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
    # Clamp shift to [0,63] before UInt64 conversion — when not subnormal,
    # shift_sub is negative but the result is overridden by the select chain.
    shift_clamped = clamp(shift_sub, Int64(0), Int64(63))
    shift_u = UInt64(ifelse(flush_to_zero, Int64(0), shift_clamped))
    lost_mask_sub = (UInt64(1) << shift_u) - UInt64(1)
    lost_sub = ifelse((wr & lost_mask_sub) != UInt64(0), UInt64(1), UInt64(0))
    wr_sub_result = (wr >> shift_u) | lost_sub
    flushed_result = result_sign << 63          # flush to signed zero

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
    mant_overflow = frac_rounded == IMPLICIT    # mantissa overflow from rounding
    frac_final = ifelse(round_up,
                 ifelse(mant_overflow, UInt64(0), frac_rounded),
                 frac)
    exp_after_round = ifelse(round_up & mant_overflow,
                             result_exp + Int64(1),
                             result_exp)
    exp_overflow_after_round = exp_after_round >= Int64(0x7FF)

    # ── Pack normal result ──
    # Clamp exponent to valid UInt64 range for packing; out-of-range cases are
    # overridden by the select chain below (subnormal/overflow/flush).
    exp_pack = UInt64(clamp(exp_after_round, Int64(0), Int64(0x7FE)))
    normal_result = (result_sign << 63) | (exp_pack << 52) | frac_final

    # ── Final select chain: priority order ──
    # NaN > Inf > Zero > Subnormal flush > Exp overflow > Normal
    result = normal_result
    result = ifelse(exp_overflow | exp_overflow_after_round, overflow_result, result)
    result = ifelse(subnormal & flush_to_zero, flushed_result, result)
    result = ifelse(exact_cancel, UInt64(0), result)
    result = ifelse(a_zero & b_zero, zero_zero_result, result)
    result = ifelse(b_zero & (!a_zero), a, result)
    result = ifelse(a_zero & (!b_zero), b, result)
    result = ifelse(a_inf & b_inf, inf_inf_result, result)
    result = ifelse(b_inf & (!a_inf), inf_finite_b, result)
    result = ifelse(a_inf & (!b_inf), inf_finite_a, result)
    result = ifelse(a_nan | b_nan, QNAN, result)

    return result
end
