#!/usr/bin/env julia
"""
Bennett.jl benchmark suite — generates BENCHMARKS.md with gate counts
compared against published results.

Usage: julia --project=. benchmark/run_benchmarks.jl
"""

using Bennett
using Bennett: simulate, verify_reversibility, gate_count, ancilla_count,
               t_count, t_depth, peak_live_wires

# ---- benchmark infrastructure ----

struct BenchResult
    name::String
    width::String
    total::Int
    not_gates::Int
    cnot_gates::Int
    toffoli_gates::Int
    wires::Int
    ancillae::Int
    t_count::Int
    published_ref::String  # "Author Year: N gates" or ""
end

results = BenchResult[]

function bench!(name, f, types...; published="", width="")
    Bennett._reset_names!()
    c = reversible_compile(f, types...)
    gc = gate_count(c)
    tc = t_count(c)
    ac = ancilla_count(c)
    w = isempty(width) ? join(["$(sizeof(T)*8)" for T in types], "×") : width
    push!(results, BenchResult(name, w, gc.total, gc.NOT, gc.CNOT, gc.Toffoli,
                               c.n_wires, ac, tc, published))
    # Verify correctness
    @assert verify_reversibility(c) "Reversibility check failed for $name"
end

# ---- integer arithmetic ----

println("Running integer benchmarks...")

for (W, T) in [(8, Int8), (16, Int16), (32, Int32), (64, Int64)]
    f_add(x) = x + one(T)
    bench!("x+1", f_add, T; width="i$W",
           published="Cuccaro 2004: $(2*W) Toff (in-place)")
end

f_poly8(x::Int8) = x * x + Int8(3) * x + Int8(1)
bench!("x²+3x+1", f_poly8, Int8; width="i8")

f_mul8(x::Int8, y::Int8) = x * y
bench!("x*y", f_mul8, Int8, Int8; width="i8×i8")

f_mul32(x::Int32, y::Int32) = x * y
bench!("x*y", f_mul32, Int32, Int32; width="i32×i32")

# Cuccaro in-place adder comparison
println("Running Cuccaro comparison...")
for (W, T) in [(8, Int8), (32, Int32), (64, Int64)]
    f_add(x) = x + one(T)
    Bennett._reset_names!()
    lr_inp = Bennett.lower(Bennett.extract_parsed_ir(f_add, Tuple{T}); use_inplace=true)
    c_inp = Bennett.bennett(lr_inp)
    gc = gate_count(c_inp)
    push!(results, BenchResult("x+1 (Cuccaro)", "i$W", gc.total, gc.NOT, gc.CNOT, gc.Toffoli,
                               c_inp.n_wires, ancilla_count(c_inp), t_count(c_inp),
                               "Cuccaro 2004: $(2*W) Toff"))
end

# Constant-folded polynomial
println("Running constant-folded benchmarks...")
Bennett._reset_names!()
lr_fold = Bennett.lower(Bennett.extract_parsed_ir(f_poly8, Tuple{Int8}); fold_constants=true)
c_fold = Bennett.bennett(lr_fold)
gc = gate_count(c_fold)
push!(results, BenchResult("x²+3x+1 (folded)", "i8", gc.total, gc.NOT, gc.CNOT, gc.Toffoli,
                           c_fold.n_wires, ancilla_count(c_fold), t_count(c_fold), ""))

# ---- SHA-256 sub-functions ----

println("Running SHA-256 benchmarks...")

ch(e::UInt32, f::UInt32, g::UInt32) = (e & f) ⊻ (~e & g)
maj(a::UInt32, b::UInt32, c::UInt32) = (a & b) ⊻ (a & c) ⊻ (b & c)
rotr(x::UInt32, n::Int) = (x >> n) | (x << (32 - n))
sigma0(a::UInt32) = rotr(a, 2) ⊻ rotr(a, 13) ⊻ rotr(a, 22)
sigma1(e::UInt32) = rotr(e, 6) ⊻ rotr(e, 11) ⊻ rotr(e, 25)

function sha256_round(a::UInt32, b::UInt32, c::UInt32, d::UInt32,
                      e::UInt32, f::UInt32, g::UInt32, h::UInt32,
                      k::UInt32, w::UInt32)
    t1 = h + sigma1(e) + ch(e, f, g) + k + w
    t2 = sigma0(a) + maj(a, b, c)
    new_e = d + t1
    new_a = t1 + t2
    return (new_a, new_e)
end

bench!("SHA-256 ch", ch, UInt32, UInt32, UInt32; width="3×i32",
       published="PRS15 Fig.15: 128 Toff")
bench!("SHA-256 maj", maj, UInt32, UInt32, UInt32; width="3×i32",
       published="PRS15 Fig.15: 128 Toff")
bench!("SHA-256 Σ₀", sigma0, UInt32; width="i32")
bench!("SHA-256 Σ₁", sigma1, UInt32; width="i32")
bench!("SHA-256 round", sha256_round, ntuple(_ -> UInt32, 10)...;
       width="10×i32",
       published="PRS15 Table II: 683 Toff (hand-opt)")

# SHA-256 with constant folding
Bennett._reset_names!()
lr_sha_fold = Bennett.lower(Bennett.extract_parsed_ir(sha256_round,
    Tuple{ntuple(_ -> UInt32, 10)...}); fold_constants=true)
c_sha_fold = Bennett.bennett(lr_sha_fold)
gc_sf = gate_count(c_sha_fold)
push!(results, BenchResult("SHA-256 round (folded)", "10×i32", gc_sf.total, gc_sf.NOT,
    gc_sf.CNOT, gc_sf.Toffoli, c_sha_fold.n_wires, ancilla_count(c_sha_fold),
    t_count(c_sha_fold), "PRS15 Table II: 683 Toff (hand-opt)"))

# ---- Float64 operations ----

println("Running Float64 benchmarks...")

bench!("soft_fadd", x -> x + 1.0, Float64; width="f64",
       published="Haener 2018: ~2000 Toff (no NaN/Inf)")
bench!("soft_fmul", (x, y) -> x * y, Float64, Float64; width="f64×f64")

# ---- optimization comparison ----

println("Running optimization comparisons...")

f_inc(x::Int8) = x + Int8(3)
Bennett._reset_names!()
lr_std = Bennett.lower(Bennett.extract_parsed_ir(f_inc, Tuple{Int8}))
Bennett._reset_names!()
lr_inp = Bennett.lower(Bennett.extract_parsed_ir(f_inc, Tuple{Int8}); use_inplace=true)

c_full = Bennett.bennett(lr_std)
c_inplace = Bennett.bennett(lr_inp)
c_eager = Bennett.value_eager_bennett(lr_inp)

println("\n=== Optimization comparison: x+3 (Int8) ===")
println("  Full Bennett:     $(c_full.n_wires) wires, peak=$(peak_live_wires(c_full))")
println("  Cuccaro in-place: $(c_inplace.n_wires) wires, peak=$(peak_live_wires(c_inplace))")
println("  Cuccaro+EAGER:    $(c_eager.n_wires) wires, peak=$(peak_live_wires(c_eager))")

# SHA-256 pebbled
Bennett._reset_names!()
parsed_sha = Bennett.extract_parsed_ir(sha256_round, Tuple{ntuple(_ -> UInt32, 10)...})
lr_sha = Bennett.lower(parsed_sha)
c_sha_full = Bennett.bennett(lr_sha)
c_sha_peb = pebbled_group_bennett(lr_sha; max_pebbles=Bennett.min_pebbles(length(lr_sha.gate_groups)))
println("\n=== SHA-256 round pebbling ===")
println("  Full Bennett: $(c_sha_full.n_wires) wires, $(ancilla_count(c_sha_full)) ancillae")
println("  Pebbled(s=$(Bennett.min_pebbles(length(lr_sha.gate_groups)))): $(c_sha_peb.n_wires) wires, $(ancilla_count(c_sha_peb)) ancillae")

# ---- generate BENCHMARKS.md ----

println("\nGenerating BENCHMARKS.md...")

open(joinpath(@__DIR__, "..", "BENCHMARKS.md"), "w") do io
    println(io, "# Bennett.jl Benchmarks")
    println(io)
    println(io, "Auto-generated by `benchmark/run_benchmarks.jl`. All circuits verified reversible.")
    println(io)
    println(io, "## Gate Counts")
    println(io)
    println(io, "| Function | Width | Total | NOT | CNOT | Toffoli | Wires | Ancillae | T-count | Published |")
    println(io, "|----------|-------|-------|-----|------|---------|-------|----------|---------|-----------|")
    for r in results
        pub = isempty(r.published_ref) ? "" : r.published_ref
        println(io, "| $(r.name) | $(r.width) | $(r.total) | $(r.not_gates) | $(r.cnot_gates) | $(r.toffoli_gates) | $(r.wires) | $(r.ancillae) | $(r.t_count) | $(pub) |")
    end
    println(io)
    println(io, "## Optimization Comparison")
    println(io)
    println(io, "### x+3 (Int8)")
    println(io, "| Strategy | Wires | Peak Live |")
    println(io, "|----------|-------|-----------|")
    println(io, "| Full Bennett | $(c_full.n_wires) | $(peak_live_wires(c_full)) |")
    println(io, "| Cuccaro in-place | $(c_inplace.n_wires) | $(peak_live_wires(c_inplace)) |")
    println(io, "| Cuccaro + EAGER | $(c_eager.n_wires) | $(peak_live_wires(c_eager)) |")
    println(io)
    println(io, "### SHA-256 Round")
    println(io, "| Strategy | Wires | Ancillae |")
    println(io, "|----------|-------|----------|")
    println(io, "| Full Bennett | $(c_sha_full.n_wires) | $(ancilla_count(c_sha_full)) |")
    s = Bennett.min_pebbles(length(lr_sha.gate_groups))
    println(io, "| Pebbled (s=$s) | $(c_sha_peb.n_wires) | $(ancilla_count(c_sha_peb)) |")
    println(io, "| PRS15 Table II (hand-opt) | 353 | — |")
    println(io, "| PRS15 Table II (Bennett) | 704 | — |")
    println(io, "| PRS15 Table II (EAGER) | 353 | — |")
end

println("Done! See BENCHMARKS.md")
