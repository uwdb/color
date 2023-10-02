
using Test
using Graphs

# This test suite aims to determine the correctness of the 'get_exact_size' function.
# We compare the results of the function with the exact sizes of known graphs.

@testset "exact symmetric graphs" begin

    @testset "1-edge graph, unlabeled" begin
        g = DataGraph(2)
        add_labeled_edge!(g, (1, 2), -1)
        query_graph = QueryGraph(2)
        add_labeled_edge!(query_graph, (1, 2), -1)
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 1
        @test predicted_size == exact_size
    end

    @testset "query larger than 1-edge graph" begin
         numVertices = 2
         g = DataGraph(numVertices)
         add_labeled_edge!(g, (1, 2), 1)
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), 1)
         add_labeled_edge!(query_graph, (2, 3), 1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 0
         @test predicted_size == exact_size
    end

     @testset "cycle graph, 1 edge query" begin
         numVertices = 1000
         g = DataGraph(numVertices)
         for i in 1:999
            add_labeled_edge!(g, (i, i+1), 1)
         end
         add_labeled_edge!(g, (1000, 1), 1)
         query_graph = QueryGraph(2)
         add_labeled_edge!(query_graph, (1, 2), 1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 1000
         @test predicted_size == exact_size
     end

    @testset "large turan graph, 2-edge query" begin
         numVertices = 12
         numPartitions = 4
         g = turan_graph(numVertices, numPartitions)
         g = DataGraph(DiGraph(g))
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 972
         @test predicted_size == exact_size
     end

    @testset "undirected cycle graph, 2 edge query" begin
        numVertices = 1000
        g = cycle_graph(numVertices)
        g = DataGraph(DiGraph(g))
        query_graph = QueryGraph(3)
        add_labeled_edge!(query_graph, (1, 2), -1)
        add_labeled_edge!(query_graph, (2, 3), -1)
        exact_size_partial = only(get_exact_size(query_graph, g; use_partial_sums=true, verbose=false))
        exact_size_no_partial = only(get_exact_size(query_graph, g; use_partial_sums=false, verbose=false))
        predicted_size = 4000
        @test predicted_size == exact_size_partial
        @test predicted_size == exact_size_no_partial
     end

    @testset "query larger than cycle graph" begin
         numVertices = 1000
         g = DataGraph(cycle_digraph(numVertices))
         query_graph = QueryGraph(2)
         add_labeled_edge!(query_graph, (1, 2), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 1000
         @test predicted_size == exact_size
    end

    @testset "Dorogovtsev-Mendes graph" begin
         numVertices = 4
         g = dorogovtsev_mendes(numVertices)
         g = DataGraph(DiGraph(g))
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 26
         @test predicted_size == exact_size
     end

    @testset "undirected binary tree graph" begin
        depth = 3
        g = binary_tree(depth)
        g = DataGraph(DiGraph(g))
        query_graph = QueryGraph(3)
        add_labeled_edge!(query_graph, (1, 2), -1)
        add_labeled_edge!(query_graph, (2, 3), -1)
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 26
        @test predicted_size == exact_size
     end

    @testset "directed binary tree graph" begin
         numVertices = 7
         g = DataGraph(numVertices)
         add_labeled_edge!(g, (1, 2), 1)
         add_labeled_edge!(g, (1, 3), 1)
         add_labeled_edge!(g, (2, 4), 1)
         add_labeled_edge!(g, (2, 5), 1)
         add_labeled_edge!(g, (3, 6), 1)
         add_labeled_edge!(g, (3, 7), 1)
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 4
         @test predicted_size == exact_size
    end

    @testset "star graph" begin
         numVertices = 4
         g = star_graph(numVertices)
         g = DataGraph(DiGraph(g))
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 12
         @test predicted_size == exact_size
     end

    @testset "disconnected graph" begin
         numVertices = 4
         g = DataGraph(numVertices)
         add_labeled_edge!(g, (1, 2), 1)
         add_labeled_edge!(g, (3, 4), 1)
         query_graph = QueryGraph(2)
         add_labeled_edge!(query_graph, (1, 2), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 2
         @test predicted_size == exact_size
    end


    @testset "directed triangle query on undirected triangle data" begin
        g = DataGraph(DiGraph(cycle_graph(3)))
        query_graph = QueryGraph(3)
        add_labeled_edge!(query_graph, (1,2), -1)
        add_labeled_edge!(query_graph, (2,3), -1)
        add_labeled_edge!(query_graph, (3,1), -1)
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 3*2
        @test predicted_size == exact_size
    end

    @testset "directed triangle query on undirected 4-cycle data" begin
        g = DataGraph(DiGraph(cycle_graph(4)))
        query_graph = QueryGraph(3)
        add_labeled_edge!(query_graph, (1,2), -1)
        add_labeled_edge!(query_graph, (2,3), -1)
        add_labeled_edge!(query_graph, (3,1), -1)
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 0
        @test predicted_size == exact_size
    end

    @testset "directed triangle query on undirected 10 clique data" begin
        g = DataGraph(DiGraph(clique_graph(10, 1)))
        for i in 2:10
            add_labeled_edge!(g, (i, i), -1) # The clique_graph constructor doesn't include self-edges, except for node 1...
        end
        query_graph = QueryGraph(3)
        add_labeled_edge!(query_graph, (1,2), -1)
        add_labeled_edge!(query_graph, (2,3), -1)
        add_labeled_edge!(query_graph, (3,1), -1)
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 10*10*10

        @test predicted_size == exact_size
    end

    @testset "1-vertex query" begin
         g = DataGraph(DiGraph(path_graph(2)))
         update_node_labels!(g, 1, [1])
         update_node_labels!(g, 2, [1])
         query_graph = QueryGraph(1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 2
         @test predicted_size == exact_size
     end
end

 @testset "asymmetric graphs" begin
     @testset "undirected asymmetric star graph" begin
         numVertices = 6
         g = Graph(numVertices)
         add_edge!(g, 1, 2)
         add_edge!(g, 1, 3)
         add_edge!(g, 1, 4)
         add_edge!(g, 1, 5)
         add_edge!(g, 5, 6)
         g = DataGraph(DiGraph(g))
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 24
         @test predicted_size == exact_size
     end

     @testset "undirected asymmetric barbell graph" begin
         numVertices = 7
         g = Graph(numVertices)
         add_edge!(g, 1, 2)
         add_edge!(g, 1, 3)
         add_edge!(g, 2, 4)
         add_edge!(g, 3, 4)
         add_edge!(g, 4, 5)
         add_edge!(g, 5, 7)
         add_edge!(g, 5, 6)
         add_edge!(g, 6, 7)
         g = DataGraph(DiGraph(g))
         query_graph = QueryGraph(3)
         add_labeled_edge!(query_graph, (1, 2), -1)
         add_labeled_edge!(query_graph, (2, 3), -1)
         exact_size = only(get_exact_size(query_graph, g; verbose=false))
         predicted_size = 38
         @test predicted_size == exact_size
     end
 end
