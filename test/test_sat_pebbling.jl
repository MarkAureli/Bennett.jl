@testset "SAT-based pebbling (Meuli 2019)" begin

    @testset "linear chain: 3 nodes, full Bennett" begin
        # A → B → C, output = C. Min pebbles = 3 (Knill: 1+ceil(log2(3)))
        adj = Dict(1 => Int[], 2 => [1], 3 => [2])
        outputs = [3]

        # Full Bennett (P=3): should find solution
        schedule = Bennett.sat_pebble(adj, outputs; max_pebbles=3, max_steps=5)
        @test !isnothing(schedule)
        @test Bennett.verify_pebble_schedule(adj, outputs, schedule)
        println("  chain(3), P=3: $(length(schedule)-1) steps")
    end

    @testset "diamond DAG: 4 nodes" begin
        # A,B inputs → C depends on A,B → D depends on C. Output = D.
        adj = Dict(1 => Int[], 2 => Int[], 3 => [1, 2], 4 => [3])
        outputs = [4]

        schedule = Bennett.sat_pebble(adj, outputs; max_pebbles=4, max_steps=7)
        @test !isnothing(schedule)
        @test Bennett.verify_pebble_schedule(adj, outputs, schedule)
        println("  diamond(4), P=4: $(length(schedule)-1) steps")
    end

    @testset "chain(5): SAT finds reduced-pebble schedule" begin
        # Chain of 5 nodes. Full Bennett: P=5. Min: P=4 (1+ceil(log2(5))=4)
        adj = Dict(i => (i == 1 ? Int[] : [i-1]) for i in 1:5)
        outputs = [5]

        # Full Bennett (P=5)
        full = Bennett.sat_pebble(adj, outputs; max_pebbles=5, max_steps=9)
        @test !isnothing(full)
        @test Bennett.verify_pebble_schedule(adj, outputs, full)

        # Reduced (P=4) — needs more steps (Knill: F(5,4)=7)
        reduced = Bennett.sat_pebble(adj, outputs; max_pebbles=4, timeout_steps=30)
        @test !isnothing(reduced)
        @test Bennett.verify_pebble_schedule(adj, outputs, reduced)
        println("  chain(5): P=5 → $(length(full)-1) steps, P=4 → $(length(reduced)-1) steps")
    end

    @testset "verify_pebble_schedule rejects bad schedules" begin
        adj = Dict(1 => Int[], 2 => [1], 3 => [2])
        outputs = [3]
        # Bad: skip pebbling prerequisite
        bad = [Set{Int}(), Set([2]), Set([3])]
        @test !Bennett.verify_pebble_schedule(adj, outputs, bad)
    end
end
