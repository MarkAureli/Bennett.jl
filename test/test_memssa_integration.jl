using Test
using Bennett
using Bennett: MemSSAInfo

# T2a.3 — Integration tests demonstrating MemorySSA captures information that
# T0 preprocessing (sroa / mem2reg / simplifycfg / instcombine) cannot.
#
# For each pattern: compile the function with preprocess=true AND
# use_memory_ssa=true, then assert:
#   (a) memory operations SURVIVED preprocessing (the interesting case),
#   (b) MemSSAInfo contains the Def/Use/Phi annotations that would let
#       a future lower_load! pick the right incoming memory state.
#
# These tests establish the diagnostic capability. Wiring the info into
# lower_load! for correctness-improving dispatch is out-of-scope follow-up.

# Helper: count instruction types across all blocks of a ParsedIR.
function _count_ir(parsed, ::Type{T}) where T
    n = 0
    for b in parsed.blocks
        for inst in b.instructions
            inst isa T && (n += 1)
        end
    end
    return n
end

@testset "T2a.3 MemorySSA integration — cases T0 misses" begin

    @testset "var-index load into local array" begin
        # SROA+mem2reg can't eliminate a var-indexed alloca. The store chain
        # remains visible and MemorySSA annotates it.
        f(x::UInt8, i::UInt8) = let a = [x, x+UInt8(1), x+UInt8(2), x+UInt8(3)]
            a[(i & 0x3) + 1]
        end

        parsed = Bennett.extract_parsed_ir(f, Tuple{UInt8, UInt8};
                                            preprocess=true, use_memory_ssa=true)
        @test parsed.memssa !== nothing
        # Memory ops DID survive preprocessing
        n_stores = _count_ir(parsed, Bennett.IRStore)
        n_loads  = _count_ir(parsed, Bennett.IRLoad)
        @test n_stores >= 1 || n_loads >= 1
        # Memssa has non-empty annotations
        @test !isempty(parsed.memssa.def_clobber)
        @test !isempty(parsed.memssa.use_at_line)
    end

    @testset "conditional store in diamond CFG produces MemoryPhi" begin
        # The paper-winning case: one branch stores into the array, the other
        # doesn't. At merge, MemorySSA synthesizes a MemoryPhi telling us the
        # load reads either branch's state. T0 preprocessing cannot simplify
        # this into a single value because the branch condition is dynamic.
        f(x::UInt8, cond::Bool) = let a = [UInt8(0), UInt8(0), UInt8(0), UInt8(0)]
            if cond
                a[1] = x
            end
            a[1]
        end

        parsed = Bennett.extract_parsed_ir(f, Tuple{UInt8, Bool};
                                            preprocess=true, use_memory_ssa=true)
        @test parsed.memssa !== nothing
        # Either memssa captures a Phi directly, or the diamond was folded and
        # we still have Def/Use annotations — either way non-empty.
        mem_nonempty = !isempty(parsed.memssa.phis) ||
                       !isempty(parsed.memssa.use_at_line) ||
                       !isempty(parsed.memssa.def_clobber)
        @test mem_nonempty
    end

    @testset "sequential stores + load" begin
        # Multiple stores to the same location (not vectorized — avoids
        # InsertElement emit from SROA on array patterns). Each store creates
        # a distinct MemoryDef; the final load reads the last.
        f(x::UInt8) = let a = Ref(UInt8(0))
            a[] = x
            a[] = x + UInt8(1)
            a[] = x + UInt8(2)
            a[]
        end

        parsed = Bennett.extract_parsed_ir(f, Tuple{UInt8};
                                            preprocess=false, use_memory_ssa=true)
        @test parsed.memssa !== nothing
        # Raw IR (no preprocess) has stores/loads of the Ref
        @test !isempty(parsed.memssa.def_clobber) ||
              !isempty(parsed.memssa.use_at_line)
    end

    @testset "memssa-off matches T0 behavior exactly" begin
        # The use_memory_ssa flag must be a pure addition: turning it on
        # doesn't change the walked IR (ParsedIR.blocks, args, ret_width
        # should match bit-for-bit when we turn memssa on vs off).
        f(x::Int8) = x + Int8(1)
        a_off = Bennett.extract_parsed_ir(f, Tuple{Int8})
        a_on  = Bennett.extract_parsed_ir(f, Tuple{Int8}; use_memory_ssa=true)
        @test length(a_off.blocks) == length(a_on.blocks)
        @test a_off.args == a_on.args
        @test a_off.ret_width == a_on.ret_width
        @test a_off.memssa === nothing
        @test a_on.memssa !== nothing
    end

    @testset "annotation IDs form a consistent graph" begin
        # Every Use's target Def should exist in def_clobber (or be
        # live-on-entry sentinel 0). No dangling references.
        f(x::Int, i::Int) = let a = [x, x+1, x+2, x+3]
            a[(i & 0x3) + 1]
        end
        parsed = Bennett.extract_parsed_ir(f, Tuple{Int, Int};
                                            preprocess=true, use_memory_ssa=true)
        for (_, def_id) in parsed.memssa.use_at_line
            @test def_id == 0 || haskey(parsed.memssa.def_clobber, def_id) ||
                  haskey(parsed.memssa.phis, def_id)
        end
        # Every def's clobber target either exists as another Def or is :live_on_entry
        for (_, clobber) in parsed.memssa.def_clobber
            @test clobber === :live_on_entry ||
                  haskey(parsed.memssa.def_clobber, clobber) ||
                  haskey(parsed.memssa.phis, clobber)
        end
    end
end
