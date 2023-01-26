include("../../Source/UnlabeledCardinalityEstimator.jl")
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!

using Test
using Graphs

@testset "exact symmetrical graphs" begin
    @testset "1-edge graph" begin
        numVertices = 2
        g = path_graph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2
        @test predicted_size == exact_size
    end

    @testset "query larger than 1-edge graph" begin
        numVertices = 2
        g = path_graph(numVertices)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2,3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2
        @test predicted_size == exact_size
    end
    
    @testset "cycle graph, 1 edge query" begin
        numVertices = 1000
        g = cycle_graph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2000
        @test predicted_size == exact_size
    end

    @testset "large turan graph, 2-edge query" begin
        numVertices = 12
        numPartitions = 4
        g = turan_graph(numVertices, numPartitions)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 972
        @test predicted_size == exact_size
    end

    @testset "cycle graph, 2 edge query" begin
        numVertices = 1000
        g = cycle_graph(numVertices)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 4000
        @test predicted_size == exact_size
    end

    @testset "query larger than cycle graph" begin
        numVertices = 1000
        g = cycle_graph(numVertices)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2000
        @test predicted_size == exact_size
    end

    @testset "Dorogovtsev-Mendes graph" begin
        numVertices = 4
        g = dorogovtsev_mendes(numVertices)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 26
        @test predicted_size == exact_size
    end

    @testset "binary tree graph" begin
        depth = 3
        g = binary_tree(depth)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 26
        @test predicted_size == exact_size
    end

    @testset "star graph" begin
        numVertices = 4
        g = star_graph(numVertices)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 12
        @test predicted_size == exact_size
    end

    # result is too large (5.36e303) to predict, can't be used to test
    # validity of get_exact_size function
    # @testset "simple path graph" begin
    #     numVertices = 1000
    #     g = cycle_graph(numVertices)
    #     query_graph = DiGraph(1000)
    #     for i in 1:999
    #         add_edge!(query_graph, (i, i+1))
    #     end
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     predicted_size = 2000
    #     @test predicted_size == exact_size
    # end
end

@testset "asymmetrical graphs" begin
    @testset "asymmetrical star graph" begin
        numVertices = 6
        g = Graph(numVertices)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 1, 5)
        add_edge!(g, 5, 6)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 24
        @test predicted_size == exact_size
    end

    @testset "asymmetrical barbell graph" begin
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
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2, 3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 38
        @test predicted_size == exact_size
    end
end