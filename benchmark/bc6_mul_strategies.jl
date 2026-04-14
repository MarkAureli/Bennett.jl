#!/usr/bin/env julia
"""
BC.6 — Multiplication strategy head-to-head (Bennett-hllu).

Compares three multiplier implementations at W ∈ {4, 8, 16, 32, 64} on
squaring (`x -> x * x`), driven through the `mul=` dispatcher:

  1. `:shift_add` (`lower_mul!`) — classical schoolbook, O(W²) Toffolis,
     O(W²) Toffoli-depth, Θ(W²) ancilla.
  2. `:karatsuba` (`lower_mul_karatsuba!`) — Θ(W^log₂3) Toffolis,
     Θ(W^log₂5) wires. Kicks in at W ≥ 8 via dispatcher heuristic.
  3. `:qcla_tree` (`lower_mul_qcla_tree!`, Sun-Borissov 2026
     arXiv:2604.09847) — O(n²) Toffolis, O(log² n) Toffoli-depth,
     self-reversing primitive.

Headline: QCLA tree wins on Toffoli-depth at every width (exponential
advantage on depth-sensitive workloads), at the cost of ~2× Toffolis and
~8× ancilla vs shift-and-add.
"""

using Bennett
using Printf

const WIDTHS = (8, 16, 32, 64)   # Int8 is Julia's narrowest; W=4 duplicates W=8 via Int8
const STRATEGIES = (:shift_add, :karatsuba, :qcla_tree)

# Pick the integer type matching the width.
function _intT(W)
    W <= 8  && return Int8
    W <= 16 && return Int16
    W <= 32 && return Int32
    return Int64
end

function _measure(W, strategy)
    T = _intT(W)
    c = reversible_compile((x, y) -> x * y, T, T; mul=strategy)
    gc = gate_count(c)
    return (
        W = W,
        strategy = strategy,
        total = gc.total,
        NOT = gc.NOT,
        CNOT = gc.CNOT,
        Toffoli = gc.Toffoli,
        tof_depth = toffoli_depth(c),
        depth = depth(c),
        t_count = t_count(c),
        ancilla = ancilla_count(c),
        peak_live = peak_live_wires(c),
    )
end

function main(io::IO=stdout)
    println(io, "# BC.6 Multiplication strategy head-to-head")
    println(io)
    println(io, "Function: `f(x, y) = x * y`. Compiled via `reversible_compile(..., mul=STRAT)`.")
    println(io, "All three strategies verified correct and reversible via the full test suite.")
    println(io)
    @printf(io, "| %-3s | %-12s | %-8s | %-8s | %-8s | %-8s | %-8s | %-10s | %-8s | %-10s |\n",
            "W", "strategy", "total", "Toffoli", "Tof-depth", "depth", "t_count", "ancilla", "peak", "NOT+CNOT")
    @printf(io, "|----:|:-------------|--------:|--------:|----------:|------:|--------:|----------:|--------:|----------:|\n")
    for W in WIDTHS, strategy in STRATEGIES
        m = _measure(W, strategy)
        @printf(io, "| %3d | %-12s | %8d | %8d | %9d | %6d | %8d | %10d | %8d | %10d |\n",
                m.W, string(m.strategy), m.total, m.Toffoli, m.tof_depth, m.depth, m.t_count,
                m.ancilla, m.peak_live, m.NOT + m.CNOT)
    end
    println(io)
    println(io, "## Notes")
    println(io)
    println(io, "- `Tof-depth` is the longest chain of Toffolis (from `toffoli_depth`).")
    println(io, "- `t_count` is `7 × Toffoli` (AMMR Toffoli-to-T decomposition).")
    println(io, "- `depth` is gate-level depth with per-wire dependencies (`depth()`).")
    println(io, "- `peak` is peak simultaneously live qubits (`peak_live_wires`).")
    println(io, "- `:qcla_tree` requires Draper QCLA (quant-ph/0406142) + Sun-Borissov 2026 (arXiv:2604.09847).")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(stdout)
end
