include("../src/QuasiStableCardinalityEstimator.jl")
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!

using Test
using Graphs

# This test suite aims to determine the correctness of the 'sample_paths' function.
# We take a list of predefined partial paths and make sure the sampled paths will appropriately
# redistribute weight among the sampled paths

@testset "simple cases" begin

    @testset "identical partial path bounds" begin
        partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}} = []
        push!(partial_paths, ([1, 1], [1, 2, 3]))
        push!(partial_paths, ([2, 2], [1, 2, 3]))
        push!(partial_paths, ([3, 3], [1, 2, 3]))
        push!(partial_paths, ([4, 4], [1, 2, 3]))

        num_samples = 2

        sampled_paths = sample_paths(partial_paths, num_samples)

        for path_and_bounds in sampled_paths
            @test path_and_bounds[2][1] == 2
            @test path_and_bounds[2][2] == 4
            @test path_and_bounds[2][3] == 6
        end
    end
end
