include("../Source/UnlabeledCardinalityEstimator.jl")
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!

using Test
using Graphs

# This test suite aims to determine the correctness of the 'get_cardinality_bounds' function.
# We compare the results of the function with the exact sizes and make sure the min/avg/max are appropriate.
# We also test to make sure partial aggregation has a minimal effect on the results.

@testset "exact symmetrical graphs" begin

    @testset "1-edge graph" begin
        g = DiGraph(2)
        add_edge!(g, (1, 2))
        g_edge_labels::Dict{Int, Dict{Int, Array{Int}}} = Dict()
        g_edge_labels[1] = Dict()
        g_edge_labels[1][2] = Array([1])
        g_vertex_labels::Dict{Int, Array{Int}} = Dict()
        g_vertex_labels[1] = Array([1])
        g_vertex_labels[2] = Array([1])
        g_property = PropertyGraph(g, g_edge_labels, g_vertex_labels)
        summary = generate_color_summary(g_property, 16)
        query_graph = DiGraph(2)
        add_edge!(query_graph, (1, 2))
        q_edge_labels::Dict{Int, Dict{Int, Array{Int}}} = Dict()
        q_edge_labels[1] = Dict()
        q_edge_labels[1][2] = Array([1])
        q_vertex_labels::Dict{Int, Array{Int}} = Dict()
        q_vertex_labels[1] = Array([1])
        q_vertex_labels[2] = Array([1])
        q_property = PropertyGraph(query_graph, q_edge_labels, q_vertex_labels)
        exact_size = only(get_exact_size(q_property, g_property; verbose=false))
        bounds_without_partial_agg = get_cardinality_bounds(q_property, summary; use_partial_sums = false, verbose = false);
        bounds_with_partial_agg = get_cardinality_bounds(q_property, summary; use_partial_sums = true, verbose = false);
        # test that min/avg/max are reasonable for bounds without partial sums
        @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
        @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
        @test bounds_without_partial_agg[1] <= exact_size
        @test exact_size <= bounds_without_partial_agg[3]
        # test that min/avg/max are reasonable for bounds with apartial sums
        @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
        @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
        @test bounds_with_partial_agg[1] <= exact_size
        @test exact_size <= bounds_with_partial_agg[3]
        # test that partial aggregation doesn't affect results
        @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
        @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
        @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    end

    # @testset "query larger than 1-edge graph" begin
    #     numVertices = 2
    #     g = path_digraph(numVertices)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(3)
    #     add_edge!(query_graph, (1, 2))
    #     add_edge!(query_graph, (2,3))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end
    
    # @testset "cycle graph" begin
    #     numVertices = 1000
    #     g = cycle_digraph(numVertices)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(2)
    #     add_edge!(query_graph, (1, 2))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "query larger than cycle graph" begin
    #     numVertices = 1000
    #     g = cycle_digraph(numVertices)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(2)
    #     add_edge!(query_graph, (1, 2))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "simple path graph" begin
    #     numVertices = 1000
    #     g = cycle_digraph(numVertices)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(1000)
    #     for i in 1:999
    #         add_edge!(query_graph, (i, i+1))
    #     end
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "Dorogovtsev-Mendes graph" begin
    #     numVertices = 4
    #     g = dorogovtsev_mendes(numVertices)
    #     g = DiGraph(g)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(2)
    #     add_edge!(query_graph, (1, 2))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "binary tree graph" begin
    #     numVertices = 7
    #     g = DiGraph(numVertices)
    #     add_edge!(g, (1, 2))
    #     add_edge!(g, (1, 3))
    #     add_edge!(g, (2, 4))
    #     add_edge!(g, (2, 5))
    #     add_edge!(g, (3, 6))
    #     add_edge!(g, (3, 7))
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(3)
    #     add_edge!(query_graph, (1, 2))
    #     add_edge!(query_graph, (2, 3))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "star graph" begin
    #     numVertices = 4
    #     g = DiGraph(numVertices)
    #     add_edge!(g, (1, 2))
    #     add_edge!(g, (1, 3))
    #     add_edge!(g, (1, 4))
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(3)
    #     add_edge!(query_graph, (1, 2))
    #     add_edge!(query_graph, (2, 3))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end

    # @testset "disconnected graph" begin
    #     numVertices = 4
    #     g = DiGraph(numVertices)
    #     add_edge!(g, 1, 2)
    #     add_edge!(g, 3, 4)
    #     summary = generate_color_summary(g, 16)
    #     query_graph = DiGraph(2)
    #     add_edge!(query_graph, (1, 2))
    #     exact_size = only(get_exact_size(query_graph, g; verbose=false))
    #     bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
    #     bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
    #     # test that min/avg/max are reasonable for bounds without partial sums
    #     @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
    #     @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
    #     @test bounds_without_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_without_partial_agg[3]
    #     # test that min/avg/max are reasonable for bounds with apartial sums
    #     @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
    #     @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
    #     @test bounds_with_partial_agg[1] <= exact_size
    #     @test exact_size <= bounds_with_partial_agg[3]
    #     # test that partial aggregation doesn't affect results
    #     @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
    #     @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
    #     @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
    # end
end

# @testset "asymmetrical graphs" begin
#     @testset "asymmetrical star graph" begin
#         numVertices = 6
#         g = DiGraph(numVertices)
#         add_edge!(g, 1, 2)
#         add_edge!(g, 1, 3)
#         add_edge!(g, 1, 4)
#         add_edge!(g, 1, 5)
#         add_edge!(g, 5, 6)
#         summary = generate_color_summary(g, 16)
#         query_graph = DiGraph(3)
#         add_edge!(query_graph, (1, 2))
#         add_edge!(query_graph, (2, 3))
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with apartial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end

#     @testset "asymmetrical barbell graph" begin
#         numVertices = 7
#         g = DiGraph(numVertices)
#         add_edge!(g, 1, 2)
#         add_edge!(g, 1, 3)
#         add_edge!(g, 2, 4)
#         add_edge!(g, 3, 4)
#         add_edge!(g, 4, 5)
#         add_edge!(g, 5, 7)
#         add_edge!(g, 5, 6)
#         add_edge!(g, 6, 7)
#         summary = generate_color_summary(g, 16)
#         query_graph = DiGraph(2)
#         add_edge!(query_graph, (1, 2))
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with apartial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end
# end

# @testset "cyclic queries" begin
#     @testset "simple cyclic data and query" begin
#         numVertices = 3
#         g = cycle_digraph(numVertices)
#         summary = generate_color_summary(g, 16)
#         query_graph = DiGraph(3)
#         add_edge!(query_graph, (1, 2))
#         add_edge!(query_graph, (2, 3))
#         add_edge!(query_graph, (3, 1))
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with partial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that lower bounds are 0
#         @test bounds_with_partial_agg[1] == 0
#         @test bounds_without_partial_agg[1] == 0
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end
# end

# @testset "randomly generated graphs" begin
#     @testset "zipfian graph, simple query" begin
#         n = 2000
#         numVertices = 1000
#         zipf = [1.0/(i^.5) for i in 1:numVertices]
#         zipf = zipf ./ sum(zipf)
#         d = DiscreteNonParametric(1:numVertices, zipf)
#         x1 = rand(d, n) .% numVertices
#         x2 = rand(d, n) .% numVertices
#         g = DiGraph(numVertices)
#         for i in range(1, length(x1))
#             add_edge!(g, x1[i], x2[i])
#         end
#         summary = generate_color_summary(g, 16)
#         query_graph = DiGraph(3)
#         add_edge!(query_graph, (1, 2))
#         add_edge!(query_graph, (2, 3))
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with apartial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end

#     @testset "random regular graph, simple query" begin
#         g = random_regular_digraph(500, 5)
#         summary = generate_color_summary(g, 16)
#         query_graph = DiGraph(3)
#         add_edge!(query_graph, (1, 2))
#         add_edge!(query_graph, (2, 3))
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with apartial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end

#     @testset "zipfian graph, larger random path query" begin
#         n = 2000
#         numVertices = 1000
#         zipf = [1.0/(i^.5) for i in 1:numVertices]
#         zipf = zipf ./ sum(zipf)
#         d = DiscreteNonParametric(1:numVertices, zipf)
#         x1 = rand(d, n) .% numVertices
#         x2 = rand(d, n) .% numVertices
#         g = DiGraph(numVertices)
#         for i in range(1, length(x1))
#             add_edge!(g, x1[i], x2[i])
#         end
#         summary = generate_color_summary(g, 16)
#         # it's hard to generate random query graphs since all the nodes in the query must be connected,
#         # which is not guaranteed by Julia random graph generators.
#         # Also as this approaches 10 it takes longer time to complete :/
#         numVertices = rand(range(3, 6))
#         query_graph = path_digraph(numVertices)
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with partial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end

#     @testset "Simple Graph, Path Query" begin
#         g = SimpleDiGraph(100, 200)
#         summary = generate_color_summary(g, 16)
#         # it's hard to generate random query graphs since all the nodes in the query must be connected,
#         # which is not guaranteed by Julia random graph generators.
#         numVertices = rand(range(2, 4))
#         query_graph = path_digraph(numVertices)
#         exact_size = only(get_exact_size(query_graph, g; verbose=false))
#         bounds_without_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = false, verbose = false);
#         bounds_with_partial_agg = get_cardinality_bounds(query_graph, summary; use_partial_sums = true, verbose = false);
#         # test that min/avg/max are reasonable for bounds without partial sums
#         @test bounds_without_partial_agg[1] <= bounds_without_partial_agg[2]
#         @test bounds_without_partial_agg[2] <= bounds_without_partial_agg[3]
#         @test bounds_without_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_without_partial_agg[3]
#         # test that min/avg/max are reasonable for bounds with partial sums
#         @test bounds_with_partial_agg[1] <= bounds_with_partial_agg[2]
#         @test bounds_with_partial_agg[2] <= bounds_with_partial_agg[3]
#         @test bounds_with_partial_agg[1] <= exact_size
#         @test exact_size <= bounds_with_partial_agg[3]
#         # test that partial aggregation doesn't affect results
#         @test abs(bounds_without_partial_agg[1] - bounds_with_partial_agg[1]) <= 1
#         @test abs(bounds_without_partial_agg[2] - bounds_with_partial_agg[2]) <= 1
#         @test abs(bounds_without_partial_agg[3] - bounds_with_partial_agg[3]) <= 1
#     end
# end