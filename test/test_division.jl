@testset "Integer division and remainder" begin

    @testset "unsigned division (udiv)" begin
        function udiv_test(a::UInt8, b::UInt8)::UInt8
            return b == UInt8(0) ? UInt8(0) : div(a, b)
        end
        c = reversible_compile(udiv_test, UInt8, UInt8)
        for a in UInt8(0):UInt8(15), b in UInt8(1):UInt8(15)
            expected = div(a, b)
            got = reinterpret(UInt8, Int8(simulate(c, (a, b))))
            @test got == expected
        end
        @test verify_reversibility(c)
    end

    @testset "unsigned remainder (urem)" begin
        function urem_test(a::UInt8, b::UInt8)::UInt8
            return b == UInt8(0) ? UInt8(0) : rem(a, b)
        end
        c = reversible_compile(urem_test, UInt8, UInt8)
        for a in UInt8(0):UInt8(15), b in UInt8(1):UInt8(15)
            expected = rem(a, b)
            got = reinterpret(UInt8, Int8(simulate(c, (a, b))))
            @test got == expected
        end
        @test verify_reversibility(c)
    end

    @testset "signed division (sdiv)" begin
        function sdiv_test(a::Int8, b::Int8)::Int8
            return b == Int8(0) ? Int8(0) : div(a, b)
        end
        c = reversible_compile(sdiv_test, Int8, Int8)
        for a in Int8(-8):Int8(7), b in [Int8(-4), Int8(-2), Int8(-1), Int8(1), Int8(2), Int8(4)]
            expected = div(a, b)
            got = Int8(simulate(c, (a, b)))
            @test got == expected
        end
        @test verify_reversibility(c)
    end
end
