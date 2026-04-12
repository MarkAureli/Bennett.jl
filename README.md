# Bennett.jl

**The Enzyme of reversible computation.** An LLVM-level compiler that transforms arbitrary pure Julia functions into reversible circuits (NOT, CNOT, Toffoli) via [Bennett's 1973 construction](https://doi.org/10.1137/0218053). First reversible compiler with full LLVM `store`/`alloca` support through four specialized, automatically-dispatched memory strategies.

```julia
using Bennett

# Any pure Julia function — no special types, no annotation
f(x::Int8) = x * x + Int8(3) * x + Int8(1)

circuit = reversible_compile(f, Int8)
simulate(circuit, Int8(5))    # => 41
verify_reversibility(circuit) # => true
gate_count(circuit)           # => 872

# Float64 via branchless soft-float (bit-exact with hardware)
g(x::Float64) = x^2 + 3.0*x + 1.0
circuit_f = reversible_compile(g, Float64)

# Compile-time constant tables → Babbush-Gidney QROM (O(L) Toffoli, W-independent)
h(x::UInt8) = let tbl = (UInt8(0x63), UInt8(0x7c), UInt8(0x77), UInt8(0x7b))
    tbl[(x & UInt8(0x3)) + 1]
end
circuit_tbl = reversible_compile(h, UInt8)   # 144 gates (MUX fallback would be ~7,500)

# Mutable arrays with dynamic indexing — handled automatically
k(x::Int8, y::Int8, i::Int8) = let a = Ref(x)
    a[] = x + y
    a[]
end
circuit_mut = reversible_compile(k, Int8, Int8, Int8)
```

## What this does

Given any pure, deterministic function `f`, Bennett.jl produces a reversible circuit `(x, 0) → (x, f(x))` using only NOT, CNOT, and Toffoli gates, with all ancillae verified zero after execution. The circuit is correct by construction.

The compiler extracts LLVM IR from plain Julia via [LLVM.jl](https://github.com/maleadt/LLVM.jl)'s C API, walks the IR as typed objects, lowers each instruction to reversible gates, and applies Bennett's construction. No operator overloading, no special types — plain Julia in, reversible circuit out.

## Why

- **Quantum control**: `when(qubit) do f(x) end` — any classical computation becomes a quantum-controlled operation (motivates integration with [Sturm.jl](https://github.com/tobiasosborne/Sturm.jl))
- **Reversible hardware synthesis**: direct compilation to Toffoli networks for adiabatic/reversible CMOS
- **Space-optimized quantum oracles**: Grover, phase estimation, QSVT all need reversible implementations of classical functions

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/tobiasosborne/Bennett.jl")
```

Requires Julia 1.10+ and LLVM.jl.

## Features

### LLVM instruction coverage (38 core opcodes, ≈90% of what a pure deterministic Julia function produces)

| Category | Instructions |
|----------|-------------|
| **Terminators** (4/11) | `ret`, `br`, `switch` (cascaded), `unreachable` |
| **Integer arithmetic** (13/13) | `add`, `sub`, `mul`, `udiv`, `sdiv`, `urem`, `srem`, `shl`, `lshr`, `ashr`, `and`, `or`, `xor` |
| **Float arithmetic** (5/5) | `fadd`, `fsub`, `fmul`, `fdiv`, `fneg` (via branchless soft-float) |
| **Comparison** | `icmp` (all 10 predicates), `fcmp` (6 predicates) |
| **Control flow** | `phi`, `select`, `freeze`, bounded loops (explicit unroll) |
| **Type conversion** (8/13) | `sext`, `zext`, `trunc`, `fptoui`, `fptosi`, `uitofp`, `sitofp`, `bitcast` |
| **Aggregates** (2/2) | `insertvalue`, `extractvalue` |
| **Memory** (4/7) | `alloca`, `load`, `store`, `getelementptr` |
| **Calls** | `call` with registered callees (gate-level inlining) |

**Missing opcodes** (filed as bd issues, all P3 — not critical path): vector ops (`extractelement`, `insertelement`, `shufflevector`), exception handling (`invoke`, `landingpad`, `catchpad`, etc.), atomic ops (`atomicrmw`, `cmpxchg`, `fence`), transcendental intrinsics (`sqrt`, `sin`/`cos`, `log`/`exp`/`pow`), pointer casts (`ptrtoint`, `inttoptr`, `addrspacecast`), `frem`, `fpext`/`fptrunc`, `va_arg`.

### Memory strategies — automatically dispatched per allocation site

`_pick_alloca_strategy(shape, idx)` picks the cheapest correct lowering for every store/load:

| Strategy | When it activates | Per-store cost | Per-load cost | Paper |
|----------|-------------------|----------------|---------------|-------|
| **Shadow** | static idx, any shape | 3·W CNOT, 0 Toffoli | W CNOT, 0 Toffoli | Enzyme-adapted (Moses & Churavy 2020) |
| **MUX EXCH** | dynamic idx, (W=8, N∈{4,8}) | 7,122 / 14,026 gates | 7,514 / 9,590 gates | — |
| **QROM** | read-only global constant table | — | 4(L-1) Toffoli + O(L·W) CNOT | Babbush-Gidney 2018 §III.C |
| **Feistel hash** | reversible bijective key hash | — | 8·W Toffoli | Luby-Rackoff 1988 |

Strategies compose: a function with static-idx stores (shadow) followed by a dynamic-idx load (MUX EXCH) works end-to-end — the MUX reads the post-shadow-mutation primal state.

### Benchmark headlines

| Benchmark | Bennett.jl | Baseline | Ratio |
|-----------|------------|----------|-------|
| QROM lookup L=4, W=8 | 56 gates | MUX tree 7,514 | **134× smaller** |
| QROM lookup L=8, W=8 | 144 gates | MUX tree 9,590 | **66× smaller** |
| Shadow store W=8 | 24 CNOT | MUX EXCH 7,122 | **297× smaller** |
| Feistel hash W=32 | 480 gates | Okasaki 3-node insert 71,000 | **148× smaller** |
| MD5 full (64 steps, extrap.) | ~48k Toffoli | ReVerC 2017 eager 27.5k | 1.75× |
| SHA-256 round | 1,632 Toffoli | PRS15 Table II hand-opt 683 | 2.4× |

QROM's Toffoli count is exactly **4(L-1) and independent of W** — matches the Babbush-Gidney paper bound. Feistel is exactly **8W Toffoli** (4 rounds × 2·half-width). Shadow memory is strictly **3W CNOT / W CNOT** with zero Toffolis from the mechanism itself.

### Space optimization

Multiple Bennett construction strategies for space-time tradeoffs:

| Strategy | SHA-256 round wires | Reduction |
|----------|--------------------|-----------|
| Full Bennett | 2,049 | baseline |
| `checkpoint_bennett` | 1,761 | 14% |
| `pebbled_group_bennett` (Meuli 2019 SAT pebbling) | tunable | up to 50% |
| Cuccaro in-place adder | 1,545 | 25% |
| `value_eager_bennett` (PRS15 EAGER) | further reduction | — |

### MemorySSA integration

Opt-in analysis for memory patterns that survive T0 preprocessing:

```julia
parsed = extract_parsed_ir(f, Tuple{UInt8, Bool}; use_memory_ssa=true)
parsed.memssa isa MemSSAInfo  # Def/Use/Phi graph extracted from LLVM
```

Captures via `print<memoryssa>` pass-output parsing. Informs future lowering decisions for conditional stores, aliased pointers, and phi-merged memory states. See `docs/memory/memssa_investigation.md` for the go/no-go analysis.

### Wider types and composability

- **Int8/16/32/64**: gate count scales linearly (2× per width doubling)
- **Float64**: full IEEE 754 via branchless soft-float (bit-exact with hardware on all of add/sub/mul/div/neg/cmp/fptosi/sitofp)
- **Tuple return**: `(new_a, new_e) = sha256_round(a, ..., w)`
- **NTuple input**: pointer parameters handled via static memory flattening
- **Ref**: scalar mutable state via shadow memory
- **Mutable arrays**: `[x, y, z]` through the T3b.3 universal dispatcher
- **Controlled circuits**: `controlled(circuit)` wraps every gate with a control bit
- **Function inlining**: `register_callee!(f)` enables gate-level inlining of any pure Julia function
- **Bounded loops**: explicit unrolling via `max_loop_iterations` kwarg

## Quick start

```julia
using Bennett

# Compile a function
circuit = reversible_compile(x -> x + Int8(1), Int8)

# Simulate
simulate(circuit, Int8(42))   # => 43

# Inspect
gate_count(circuit)           # => 100
ancilla_count(circuit)        # => 76
t_count(circuit)              # => 196 (Toffoli × 7)
verify_reversibility(circuit) # => true

# Controlled version (for quantum control)
cc = controlled(circuit)
simulate(cc, true, Int8(42))  # => 43 (control = true)
simulate(cc, false, Int8(42)) # => 42 (control = false)
```

## Build & test

```bash
cd Bennett.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite runs ~10,000 assertions across ~60 test files in about 90 seconds.

## Architecture

```
Julia function          LLVM IR                Parsed IR             Reversible Circuit
──────────────         ─────────              ──────────            ──────────────────
 f(x::Int8)    ──►  code_llvm() ──► extract_parsed_ir() ──► lower() ──► bennett()
                    (LLVM.jl C API)      │                    │              │
                                         │              strategy              │
                                     preprocess         dispatch              ▼
                                     sroa/mem2reg        ▼                simulate()
                                     simplifycfg    ┌──────────┐          verify_reversibility()
                                     instcombine    │ Shadow   │          gate_count() / t_count()
                                                    │ MUX EXCH │          controlled()
                                      (optional)    │ QROM     │
                                     MemorySSA      │ Feistel  │
                                     annotations    └──────────┘
```

1. **Extract** — `extract_parsed_ir(f, arg_types)` uses LLVM.jl C API to walk IR as typed objects. Optional `preprocess=true` runs sroa/mem2reg/simplifycfg/instcombine. Optional `use_memory_ssa=true` captures MemorySSA via printer-pass stderr.
2. **Lower** — `lower(parsed_ir)` maps each instruction to reversible gates. For memory ops, `_pick_alloca_strategy` picks shadow / MUX EXCH / QROM / Feistel per op.
3. **Bennett** — `bennett(lr)` applies forward + CNOT-copy + reverse (all ancillae return to zero). Alternative: `pebbled_group_bennett` (SAT-pebbling), `value_eager_bennett` (PRS15 EAGER), `checkpoint_bennett`.
4. **Simulate** — `simulate(circuit, input)` runs bit-vector simulation with ancilla-zero verification.

## Documentation

- **[Tutorial](docs/src/tutorial.md)** — compile your first reversible circuit in 10 minutes
- **[API Reference](docs/src/api.md)** — every exported function with examples
- **[Architecture Guide](docs/src/architecture.md)** — how the compiler works internally
- **[Vision PRD](Bennett-VISION-PRD.md)** — the full v1.0 roadmap and Enzyme analogy
- **[Memory design docs](docs/memory/)** — `memssa_investigation.md`, `shadow_design.md`
- **[BENCHMARKS.md](BENCHMARKS.md)** — auto-generated head-to-head tables vs published compilers
- **[WORKLOG.md](WORKLOG.md)** — the full development log; per-task gate counts, gotchas, design decisions

## Key references

| Tag | Paper | Key result |
|-----|-------|------------|
| Bennett 1989 | Time/Space Trade-Offs for Reversible Computation | O(T^{1+ε}) time, O(S·log T) space |
| Knill 1995 | Analysis of Bennett's Pebble Game | Exact pebbling recursion |
| Cuccaro 2004 | A New Quantum Ripple-Carry Addition Circuit | In-place adder: 1 ancilla, 2n-1 Toffoli |
| Luby-Rackoff 1988 | How to Construct Pseudorandom Permutations | 4-round Feistel = bijective permutation |
| PRS15 | Parent/Roetteler/Svore, Reversible Circuit Compilation | EAGER cleanup: 5.3× reduction on SHA-2 |
| Babbush-Gidney 2018 | Encoding Electronic Spectra with Linear T Complexity | QROM: 4L-4 T, independent of W |
| Enzyme 2020 | Moses & Churavy, Automatically Synthesize Fast Gradients | LLVM-level AD — inspiration for shadow memory |
| Meuli 2019 | SAT-based Pebbling | 52.77% ancilla reduction |
| ReVerC 2017 | Parent/Roetteler/Svore CAV | Prior state-of-art reversible compiler |
| Reqomp 2024 | Paradis et al., Space-constrained Uncomputation | Up to 96% reduction via lifetime-guided |

All papers downloaded to `docs/literature/` with claims verified against text.

## Project status

Memory plan critical path **complete**: four memory strategies (MUX EXCH, QROM, Feistel, Shadow) plus universal dispatcher, plus MemorySSA ingest, plus comprehensive benchmarks. See `WORKLOG.md` session log for per-task detail.

Next focus areas (per `Bennett-VISION-PRD.md`): Sturm.jl integration for quantum control (`when(qubit) do f(x) end`), full SHA-256 benchmark (BC.3), Julia EscapeAnalysis integration (T0.3), `@linear` macro for in-place-linear functions (T2b).

## License

[AGPL-3.0](LICENSE)
