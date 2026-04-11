"""Shift-and-add multiplier: result = a * b  (mod 2^W)."""
function lower_mul!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                    a::Vector{Int}, b::Vector{Int}, W::Int)
    return lower_mul_wide!(gates, wa, a, b, W, W)
end

"""
Shift-and-add widening multiplier: result = a * b with `result_width` bits.

When result_width == W: standard mod 2^W multiplication.
When result_width == 2W: full product without truncation.
"""
function lower_mul_wide!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                         a::Vector{Int}, b::Vector{Int}, W::Int, result_width::Int)
    accum = allocate!(wa, result_width)
    for i in 1:W
        shift = i - 1
        pp = allocate!(wa, result_width)
        for k in 1:W
            dest = k + shift
            dest > result_width && break
            push!(gates, ToffoliGate(a[k], b[i], pp[dest]))
        end
        new_accum = lower_add!(gates, wa, accum, pp, result_width)
        accum = new_accum
    end
    return accum
end

"""
Karatsuba multiplier: result = a * b  (mod 2^W).

Uses recursive Karatsuba with widening sub-products.
3 multiplications of ~W/2 bits instead of W^2 schoolbook gates.
"""
function lower_mul_karatsuba!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                              a::Vector{Int}, b::Vector{Int}, W::Int)
    full = _karatsuba_wide!(gates, wa, a, b, W)
    # Truncate to W bits
    return full[1:W]
end

"""
Karatsuba widening multiply: returns 2W-bit full product of two W-bit numbers.
Recurses; base case at W ≤ 4 uses schoolbook widening.
"""
function _karatsuba_wide!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                          a::Vector{Int}, b::Vector{Int}, W::Int)
    if W <= 4
        return lower_mul_wide!(gates, wa, a, b, W, 2*W)
    end

    h = W ÷ 2
    hi_w = W - h
    out_w = 2 * W  # full product width

    a_lo = a[1:h]
    a_hi = a[h+1:W]
    b_lo = b[1:h]
    b_hi = b[h+1:W]

    # z0 = a_lo * b_lo (2h-bit full product) — recursive
    z0 = _karatsuba_wide!(gates, wa, a_lo, b_lo, h)  # 2h bits

    # z2 = a_hi * b_hi (2*hi_w-bit full product) — recursive
    z2 = _karatsuba_wide!(gates, wa, a_hi, b_hi, hi_w)  # 2*hi_w bits

    # Cross sums: (h+1)-bit each
    cross_w = hi_w + 1  # max(h, hi_w) + 1 to avoid overflow
    a_cross = allocate!(wa, cross_w)
    for i in 1:h; push!(gates, CNOTGate(a_lo[i], a_cross[i])); end
    a_hi_pad = allocate!(wa, cross_w)
    for i in 1:hi_w; push!(gates, CNOTGate(a_hi[i], a_hi_pad[i])); end
    a_sum = lower_add!(gates, wa, a_cross, a_hi_pad, cross_w)

    b_cross = allocate!(wa, cross_w)
    for i in 1:h; push!(gates, CNOTGate(b_lo[i], b_cross[i])); end
    b_hi_pad = allocate!(wa, cross_w)
    for i in 1:hi_w; push!(gates, CNOTGate(b_hi[i], b_hi_pad[i])); end
    b_sum = lower_add!(gates, wa, b_cross, b_hi_pad, cross_w)

    # z1_full = a_sum * b_sum (2*cross_w-bit full product) — recursive
    z1_full = _karatsuba_wide!(gates, wa, a_sum, b_sum, cross_w)
    prod_w = 2 * cross_w

    # z1 = z1_full - z0 - z2 (at prod_w bits)
    z0_ext = allocate!(wa, prod_w)
    for i in 1:min(2*h, prod_w); push!(gates, CNOTGate(z0[i], z0_ext[i])); end
    z1_sub1 = lower_sub!(gates, wa, z1_full, z0_ext, prod_w)

    z2_ext = allocate!(wa, prod_w)
    for i in 1:min(2*hi_w, prod_w); push!(gates, CNOTGate(z2[i], z2_ext[i])); end
    z1 = lower_sub!(gates, wa, z1_sub1, z2_ext, prod_w)

    # Assemble full product: result[0:2W-1] = z0 + (z1 << h) + (z2 << 2h)
    result = allocate!(wa, out_w)
    # z0 contributes to bits 0..2h-1
    for i in 1:min(2*h, out_w); push!(gates, CNOTGate(z0[i], result[i])); end
    # z1 << h
    z1_shifted = allocate!(wa, out_w)
    for i in 1:min(prod_w, out_w - h)
        push!(gates, CNOTGate(z1[i], z1_shifted[h + i]))
    end
    partial = lower_add!(gates, wa, result, z1_shifted, out_w)
    # z2 << 2h
    z2_shifted = allocate!(wa, out_w)
    for i in 1:min(2*hi_w, out_w - 2*h)
        push!(gates, CNOTGate(z2[i], z2_shifted[2*h + i]))
    end
    final_result = lower_add!(gates, wa, partial, z2_shifted, out_w)

    return final_result
end
