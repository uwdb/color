include("../../Source/UnlabeledCardinalityEstimator.jl")
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!

using Test
using Graphs

# This test suite aims to determine the correctness of the 'get_exact_size' function.
# We compare the results of the function with the exact sizes of known graphs.

@testset "exact symmetrical graphs" begin

    @testset "1-edge graph" begin
        numVertices = 2
        g = DiGraph(2)
        add_edge!(g, (1, 2))
        #g = path_graph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 1
        @test predicted_size == exact_size
    end

    @testset "query larger than 1-edge graph" begin
        numVertices = 2
        g = DiGraph(numVertices)
        add_edge!(g, (1, 2))
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2,3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 0
        @test predicted_size == exact_size
    end
    
    @testset "cycle graph, 1 edge query" begin
        numVertices = 1000
        g = cycle_digraph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 1000
        @test predicted_size == exact_size
    end

    @testset "large turan graph, 2-edge query" begin
        numVertices = 12
        numPartitions = 4
        g = turan_graph(numVertices, numPartitions)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 972
        @test predicted_size == exact_size
    end

    @testset "undirected cycle graph, 2 edge query" begin
        numVertices = 1000
        g = cycle_graph(numVertices)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size_partial = only(get_exact_size(query_graph, g; use_partial_sums=true, verbose=false))
        exact_size_no_partial = only(get_exact_size(query_graph, g; use_partial_sums=false, verbose=false))
        predicted_size = 4000
        @test predicted_size == exact_size_partial
        @test predicted_size == exact_size_no_partial
    end

    @testset "query larger than cycle graph" begin
        numVertices = 1000
        g = cycle_digraph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 1000
        @test predicted_size == exact_size
    end

    @testset "Dorogovtsev-Mendes graph" begin
        numVertices = 4
        g = dorogovtsev_mendes(numVertices)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 26
        @test predicted_size == exact_size
    end

    @testset "undirected binary tree graph" begin
        depth = 3
        g = binary_tree(depth)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 26
        @test predicted_size == exact_size
    end

    @testset "directed binary tree graph" begin
        numVertices = 7
        g = DiGraph(numVertices)
        add_edge!(g, (1, 2))
        add_edge!(g, (1, 3))
        add_edge!(g, (2, 4))
        add_edge!(g, (2, 5))
        add_edge!(g, (3, 6))
        add_edge!(g, (3, 7))
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 4
        @test predicted_size == exact_size
    end

    @testset "star graph" begin
        numVertices = 4
        g = star_graph(numVertices)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 12
        @test predicted_size == exact_size
    end

    @testset "disconnected graph" begin
        numVertices = 4
        g = DiGraph(numVertices)
        add_edge!(g, 1, 2)
        add_edge!(g, 3, 4)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2
        @test predicted_size == exact_size
    end

    # TODO: test cyclic queries

    # The current implementation does not support 1-vertex queries, but this
    # will be fixed in the update incorporating edge/vertex labels
    # @testset "1-vertex query" begin
    #     g = path_graph(2)
    #     query_graph = DiGraph(1)
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     predicted_size = 2
    #     @test predicted_size == exact_size
    # end
end

@testset "asymmetrical graphs" begin
    @testset "undirected asymmetrical star graph" begin
        numVertices = 6
        g = Graph(numVertices)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 1, 5)
        add_edge!(g, 5, 6)
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 24
        @test predicted_size == exact_size
    end

    @testset "undirected asymmetrical barbell graph" begin
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
        g = DiGraph(g)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 38
        @test predicted_size == exact_size
    end
end