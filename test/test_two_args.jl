@testset "Two args: m(x::Int8, y::Int8) = x*y + x - y" begin
    m(x::Int8, y::Int8) = x * y + x - y
    circuit = reversible_compile(m, Int8, Int8)

    for x in Int8(0):Int8(15), y in Int8(0):Int8(15)
        @test simulate(circuit, (x, y)) == m(x, y)
    end

    @test verify_reversibility(circuit)
    println("  Two args: ", gate_count(circuit))
    print_circuit(circuit)
end
