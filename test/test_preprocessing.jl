using Test
using Bennett

@testset "T0.1 LLVM pass pipeline control in extract_parsed_ir" begin

    @testset "backward compatibility: default passes=nothing" begin
        f(x::Int8) = x + Int8(1)
        parsed_default = Bennett.extract_parsed_ir(f, Tuple{Int8})
        parsed_nothing = Bennett.extract_parsed_ir(f, Tuple{Int8}; passes=nothing)
        @test length(parsed_default.blocks) == length(parsed_nothing.blocks)
        for (b1, b2) in zip(parsed_default.blocks, parsed_nothing.blocks)
            @test length(b1.instructions) == length(b2.instructions)
        end
    end

    @testset "custom passes run without error" begin
        f(x::Int8) = x * Int8(3) + Int8(1)
        parsed = Bennett.extract_parsed_ir(f, Tuple{Int8}; passes=["sroa", "mem2reg"])
        @test parsed isa Bennett.ParsedIR
        @test !isempty(parsed.blocks)
    end

    @testset "combine with optimize=false" begin
        f(x::Int8) = x + Int8(2)
        parsed = Bennett.extract_parsed_ir(f, Tuple{Int8}; optimize=false, passes=["mem2reg"])
        @test parsed isa Bennett.ParsedIR
    end

    @testset "empty pass list is a no-op" begin
        f(x::Int8) = x + Int8(1)
        parsed_empty = Bennett.extract_parsed_ir(f, Tuple{Int8}; passes=String[])
        parsed_default = Bennett.extract_parsed_ir(f, Tuple{Int8})
        @test length(parsed_empty.blocks) == length(parsed_default.blocks)
    end

    @testset "reversible_compile still works (full pipeline)" begin
        # Backward compat: the top-level reversible_compile path should be unaffected.
        c = reversible_compile(x -> x + Int8(3), Int8)
        for x in Int8(-128):Int8(127)
            @test simulate(c, x) == x + Int8(3)
        end
        @test verify_reversibility(c)
    end
end
