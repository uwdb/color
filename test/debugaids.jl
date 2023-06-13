using Graphs: DiGraph
using StatsPlots
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!
using Test
using Graphs
include("../Source/datasets.jl")
include("../Source/QuasiStableCardinalityEstimator.jl")

# separate file so I can debug through the results of the aids estimation
aids_data_file_path = "dataset/aids/aids.txt"
        human_data_file_path = "dataset/human/human.txt"
        lubm80_data_file_path = "dataset/lubm80/lubm80.txt"
        yago_data_file_path = "dataset/yago/yago.txt"

        human_data = load_dataset(human_data_file_path)
        aids_data = load_dataset(aids_data_file_path)   

        datasets = Dict("aids"=>aids_data, "human"=>human_data)
        dataset_names = ["aids"]
        num_sample_nodes = Dict("aids"=>20, "human"=>100)
        max_cycle_size = Dict("aids"=>6, "human"=>3)

        build_time = Dict()
        summary_size = Dict()
        color_summaries = Dict()
        for dataset in dataset_names
            results = @timed generate_color_summary(datasets[dataset], 32, verbose=true, max_size = max_cycle_size[dataset], num_sample_nodes=num_sample_nodes[dataset])
            build_time[dataset] = results[2]
            summary_size[dataset] = get_color_summary_size(results[1])
            color_summaries[dataset] = results[1]
        end

        aids_query_directories = [
        "/queryset/aids/Cycle_3/",
        "/queryset/aids/Cycle_6/",
        "/queryset/aids/Flower_6/",
        "/queryset/aids/Flower_9/",
        "/queryset/aids/Flower_12/",
        "/queryset/aids/Petal_6/",
        "/queryset/aids/Petal_9/",
        "/queryset/aids/Petal_12/"]
        aids_query_paths = [readdir(pwd() * dir, join=true) for dir in aids_query_directories]
        aids_query_paths = [(aids_query_paths...)...]
        aids_exact_sizes = []
        aids_bounds = []
        aids_bounds_with_stats = []
        aids_relative_errors = []
        aids_relative_errors_with_stats = []
        aids_query_types = []
        println("Summary Size: ", summary_size["aids"])
        println("Summary Build Time: ", build_time["aids"])
        for query_path in aids_query_paths
            #println("Query: ", query_path)
            id_and_query = load_query(query_path)
            id = id_and_query[1]
            query = id_and_query[2]
            query_type = match(r".*/queryset/aids/(.*)_.*/.*", query_path).captures[1]
            bound_results = @timed get_cardinality_bounds(query, color_summaries["aids"], usingStoredStats=false)
            bound_results_with_stats = @timed get_cardinality_bounds(query, color_summaries["aids"], usingStoredStats=true)
            gcare_size = load_true_cardinality(replace(query_path, "queryset"=>"TrueCardinalities"))
            bound_results[1][2] = max(1, bound_results[1][2])
            bound_results_with_stats[1][2] = max(1, bound_results_with_stats[1][2])
            push!(aids_exact_sizes, gcare_size)
            push!(aids_bounds, bound_results[1])
            push!(aids_bounds_with_stats, bound_results_with_stats[1])
            push!(aids_relative_errors, bound_results[1] ./ gcare_size)
            push!(aids_relative_errors_with_stats, bound_results_with_stats[1] ./ gcare_size)
            push!(aids_query_types, query_type)
        end