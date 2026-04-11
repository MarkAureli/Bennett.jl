@testset "Dependency DAG extraction" begin

    @testset "simple addition DAG" begin
        f(x::Int8) = x + Int8(1)
        c = reversible_compile(f, Int8)
        dag = Bennett.extract_dep_dag(c)

        # DAG should have nodes for each gate
        @test length(dag.nodes) > 0

        # Each node should have predecessor list
        for node in dag.nodes
            @test isa(node.preds, Vector{Int})
        end

        # Output nodes should be identifiable
        @test !isempty(dag.output_nodes)
    end

    @testset "polynomial DAG" begin
        g(x::Int8) = x * x + Int8(3) * x + Int8(1)
        c = reversible_compile(g, Int8)
        dag = Bennett.extract_dep_dag(c)

        # More complex function = more nodes
        @test length(dag.nodes) > 10
        println("  polynomial DAG: $(length(dag.nodes)) nodes")
    end
end
