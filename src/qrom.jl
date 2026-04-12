# ---- QROM: Quantum Read-Only Memory for compile-time-constant tables ----
#
# Babbush-Gidney 2018 (arXiv:1805.03662v2), Section III.C Figure 10.
#
# Transformation: (idx, 0^W) → (idx, data[idx])
# where `data` is a compile-time-constant array of L words.
#
# Cost: 2(L-1) Toffoli + O(L·W) CNOT, T-count 4(L-1) independent of W.
#
# Construction (unary iteration, §III.A Fig 7): a complete binary tree of AND
# operations over log₂(L) index bits produces L leaf-flags; exactly one leaf-flag
# is active (= 1) at runtime, corresponding to the current idx. Data is encoded
# as data-dependent CNOTs from each leaf-flag into the output register. After the
# CNOT fan-out, the AND tree is reversed to uncompute all flags — self-clean.
#
# Unlike the MUX-tree used by `soft_mux_load_*` (O(L·W) gates per query — one
# ifelse per output bit per index level), QROM's T cost is independent of W.
# For a 64-bit S-box lookup at L=16, MUX is ~50k gates; QROM is ~80.

"""
    emit_qrom!(gates, wa, data::Vector{UInt64}, idx_wires, W::Int) -> Vector{Int}

Append gates to `gates` implementing the transformation `(idx, 0^W) → (idx, data[idx])`,
where `idx` is held in `idx_wires` (log₂(L) wires, LSB first) and `data` is a
compile-time-constant table of L words (W-bit each, zero-padded in UInt64).

Allocates the W output wires from `wa` and returns them. Uses log₂(L) additional
ancilla wires for the unary-iteration tree; all ancillae return to zero after
this circuit (self-uncomputing — a compute-uncompute pair of AND trees straddling
the data fan-out).

Requirements:
  * `L = length(data)` must be a power of two (≥ 1)
  * `length(idx_wires) == log₂(L)`
  * `1 ≤ W ≤ 64`

References:
  * Babbush, Gidney, Berry, Wiebe, McClean, Paler, Fowler, Neven (2018),
    "Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity",
    §III.A (unary iteration), §III.C (QROM). arXiv:1805.03662v2.
"""
function emit_qrom!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                    data::Vector{UInt64}, idx_wires::Vector{Int}, W::Int)
    L = length(data)
    L >= 1 || error("emit_qrom!: data must have ≥ 1 entry")
    1 <= W <= 64 || error("emit_qrom!: W must be in 1..64, got $W")
    n = L == 1 ? 0 : (Int(log2(L)))
    L == 1 << n || error("emit_qrom!: L must be a power of two, got L=$L")
    length(idx_wires) == n || error("emit_qrom!: idx_wires must have $n wires for L=$L, got $(length(idx_wires))")

    # Bit-width guard: data words must fit in W bits
    mask = W == 64 ? typemax(UInt64) : (UInt64(1) << W) - UInt64(1)
    for (i, d) in enumerate(data)
        (d & ~mask) == 0 || error("emit_qrom!: data[$i] = $(repr(d)) exceeds W=$W bits")
    end

    data_out = allocate!(wa, W)

    if L == 1
        # Degenerate: no index, just CNOT the single word directly onto output.
        # No flag needed — the single data word is always selected.
        word = data[1]
        for bit in 0:W-1
            if (word >> bit) & 1 == 1
                # Constant-1 output bit: apply NOT to the (zero) output wire.
                push!(gates, NOTGate(data_out[bit + 1]))
            end
        end
        return data_out
    end

    # Allocate the root control flag and initialize it to 1 (the root is always active).
    root_flag = allocate!(wa, 1)[1]
    push!(gates, NOTGate(root_flag))

    # Recursive unary iteration: walk the binary tree depth-first. At each
    # internal node we allocate two child flags (left/right), compute them
    # from the parent flag and the index bit, recurse, then uncompute.
    _qrom_tree!(gates, wa, data, data_out, idx_wires, root_flag, 0, L, n, W)

    # Uncompute root flag (reverse the initial NOT).
    push!(gates, NOTGate(root_flag))

    return data_out
end

# Depth-first walk of the binary tree. `parent_flag` is 1 iff idx ∈ [lo, hi).
# `bit_level` is the idx bit used to split this range (n-1 at root, 0 at leaves).
function _qrom_tree!(gates::Vector{ReversibleGate}, wa::WireAllocator,
                     data::Vector{UInt64}, data_out::Vector{Int},
                     idx_wires::Vector{Int}, parent_flag::Int,
                     lo::Int, hi::Int, bit_level::Int, W::Int)
    if hi - lo == 1
        # Leaf: fan parent_flag out to data_out via CNOTs per bit of data[lo+1]
        word = data[lo + 1]
        @inbounds for bit in 0:W-1
            if (word >> bit) & 1 == 1
                push!(gates, CNOTGate(parent_flag, data_out[bit + 1]))
            end
        end
        return
    end

    # Internal node: split [lo, hi) at midpoint using idx bit (bit_level - 1).
    mid = (lo + hi) >> 1
    idx_bit_wire = idx_wires[bit_level]  # idx bit (bit_level-1) is at 1-based index bit_level

    # right_flag = parent_flag AND idx_bit  (1 Toffoli)
    right_flag = allocate!(wa, 1)[1]
    push!(gates, ToffoliGate(parent_flag, idx_bit_wire, right_flag))

    # left_flag = parent_flag AND (NOT idx_bit) = parent_flag XOR right_flag
    # Implemented as: left_flag ← 0; CNOT(parent, left); CNOT(right, left).
    left_flag = allocate!(wa, 1)[1]
    push!(gates, CNOTGate(parent_flag, left_flag))
    push!(gates, CNOTGate(right_flag, left_flag))

    # Recurse — left subtree first, then right (DFS; order doesn't affect correctness).
    _qrom_tree!(gates, wa, data, data_out, idx_wires, left_flag, lo, mid, bit_level - 1, W)
    _qrom_tree!(gates, wa, data, data_out, idx_wires, right_flag, mid, hi, bit_level - 1, W)

    # Uncompute flags (reverse the three compute gates, self-inverse).
    push!(gates, CNOTGate(right_flag, left_flag))
    push!(gates, CNOTGate(parent_flag, left_flag))
    push!(gates, ToffoliGate(parent_flag, idx_bit_wire, right_flag))

    # After uncomputation both ancilla wires are zero; return them to the pool
    # so subsequent allocations reuse them (keeps wire count at O(log L) peak).
    free!(wa, [left_flag, right_flag])
    return
end
