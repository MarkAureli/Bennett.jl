using Test
using Bennett

@testset "Narrow bit-width compilation" begin
    @testset "Int4 (4-bit) addition" begin
        c = reversible_compile(x -> x + Int8(1), Int8; bit_width=4)

        # Exhaustive: all 16 inputs
        for x in 0:15
            expected = mod(x + 1, 16)
            @test simulate(c, Int8(x)) == expected
        end
        @test verify_reversibility(c)

        gc = gate_count(c)
        println("  Int4 x+1: $(gc.total) gates, $(c.n_wires) wires")
        @test gc.total < 60  # should be ~48, well under Int8's 100
    end

    @testset "Int2 (2-bit) exhaustive" begin
        c = reversible_compile(x -> x + Int8(1), Int8; bit_width=2)

        @test simulate(c, Int8(0)) == 1
        @test simulate(c, Int8(1)) == 2
        @test simulate(c, Int8(2)) == 3
        @test simulate(c, Int8(3)) == 0  # wraps mod 4
        @test verify_reversibility(c)

        gc = gate_count(c)
        println("  Int2 x+1: $(gc.total) gates, $(c.n_wires) wires")
        @test gc.total < 30  # should be ~22
    end

    @testset "Int4 polynomial" begin
        c = reversible_compile(x -> x * x + Int8(3) * x + Int8(1), Int8; bit_width=4)

        f(x) = mod(x * x + 3 * x + 1, 16)
        for x in 0:15
            @test simulate(c, Int8(x)) == f(x)
        end
        @test verify_reversibility(c)

        gc = gate_count(c)
        println("  Int4 poly: $(gc.total) gates, $(c.n_wires) wires")
        @test gc.total < gate_count(reversible_compile(x -> x * x + Int8(3) * x + Int8(1), Int8)).total
    end

    @testset "Int4 two-arg addition" begin
        c = reversible_compile((x, y) -> x + y, Int8, Int8; bit_width=4)

        for x in 0:3, y in 0:3
            @test simulate(c, (Int8(x), Int8(y))) == mod(x + y, 16)
        end
        @test verify_reversibility(c)
        println("  Int4 x+y: $(gate_count(c).total) gates")
    end

    @testset "Gate count scaling with bit_width" begin
        for W in [2, 3, 4, 6, 8]
            c = reversible_compile(x -> x + Int8(1), Int8; bit_width=W)
            gc = gate_count(c)
            println("  Int$W x+1: $(gc.total) gates")
        end
    end
end
