using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!

using Test
using Graphs

@testset "symmetrical graphs" begin
    @testset "1-edge graph" begin
        numVertices = 2
        g = path_graph(numVertices)
        summary = generate_color_summary(g, 16)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2
        @test predicted_size == exact_size
    end

    @testset "query larger than 1-edge graph" begin
        numVertices = 2
        g = path_graph(numVertices)
        summary = generate_color_summary(g, 16)
        query_graph = DiGraph(3)
        add_edge!(query_graph, (1, 2))
        add_edge!(query_graph, (2,3))
        exact_size = only(get_exact_size(query_graph, g; verbose=false))
        predicted_size = 2
        @test predicted_size == exact_size
    end
    
    @testset "cycle graph" begin
        
    end

    @testset "query larger than cycle graph" begin
        
    end

    @testset "simple path graph" begin
        
    end

    @testset "Dorogovtsev-Mendes graph" begin
        
    end

    @testset "barbell graph" begin
        
    end

    @testset "star graph" begin
        
    end
end

@testset "asymmetrical graphs" begin
    @testset "asymmetrical star graph" begin
        
    end

    @testset "asymmetrical barbell graph" begin
        
    end
end