using Plots.PlotMeasures
using Graphs
using Random
include("../Experiments.jl")

# The goal of this file is to determine the performance of just the edge update logic.
# We create a summary using only some of the edges (but all of the nodes) in the graph and then update the summary
# with the remaining edges, repeating for different proportions of the graph to update with.
# Then, we can check metrics like relative error and see how it performs with more updates.

datasets::Vector{DATASET} = [yeast, dblp]
proportions_updated = [0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, proportion_updated=current_proportion)
                                                    for current_dataset in datasets for current_proportion in proportions_updated]
println("started building")
shuffled_edges = Dict()
for experiment_params in experiment_params_list
    build_times = [("Dataset", "Partitioner", "NumColors", "BuildPhase", "BuildTime", "MemoryFootprint")]
    dataset = experiment_params.dataset
    summary_params = experiment_params.summary_params
    data = load_dataset(dataset)
    cloned_data = DataGraph(nv(data.graph))
    cloned_data.vertex_labels = data.vertex_labels
    edges_for_later = []
    shuffled_edges[dataset] = shuffle!(collect(edges(data.graph)))
    edges_to_add = ne(data.graph) - round(convert(Float64, ne(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_updated))
    for edge in shuffled_edges[dataset]
        if edges_to_add > 0
            edges_to_add -= 1
            add_labeled_edge!(cloned_data, (src(edge), dst(edge)), only(data.edge_labels[(src(edge), dst(edge))]))
        else
            push!(edges_for_later, edge)
        end
    end
    summary_name = params_to_summary_filename(experiment_params)
    summary_file_location = "Experiments/SerializedSummaries/" * summary_name
    println("Building Color Summary: ", summary_name)
    timing_vec = Float64[]
    results = @timed generate_color_summary((experiment_params.summary_params.proportion_updated > 0) ? cloned_data : data, summary_params; verbose=1, timing_vec=timing_vec)
    current_summary = results.value
    if (experiment_params.summary_params.proportion_updated > 0)
        for edge in edges_for_later
            add_summary_edge!(current_summary, src(edge), dst(edge), get(data.edge_labels, (src(edge), dst(edge)), []))
        end
    end
    summary_size = Base.summarysize(current_summary)
    serialize(summary_file_location, current_summary)
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "FullTime",
        string(results.time),
        string(summary_size)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "Coloring",
        string(timing_vec[1]),
        string(summary_size)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "CycleCounting",
        string(timing_vec[2]),
        string(summary_size)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "BloomFilter",
        string(timing_vec[3]),
        string(summary_size)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "CardinalityCounting",
        string(timing_vec[4]),
        string(summary_size)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "EdgeStats",
        string(timing_vec[5]),
        string(summary_size)))
    results_filename = params_to_results_filename(experiment_params)
    result_file_location = "Experiments/Results/Build_" * results_filename
    writedlm(result_file_location, build_times, ",")
end

println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=build_time, ylims=[0, 10], y_ticks = [0, 2, 4 ,6 ,8, 10], legend_pos=:topright, x_label="Proportion Updated", y_label="Build Time (S)", grouping=proportion_updated, filename="just-edge-updates-build")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error,ylims=[10^-20, 10^15],  x_label="Proportion Updated", y_label="Estimate Error", grouping=proportion_updated, filename="just-edge-updates-error")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=runtime, ylims=[10^-5, 10], y_ticks = [10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1, 10], x_label="Proportion Updated", y_label="Runtime (S)", grouping=proportion_updated, filename="just-edge-updates-runtime")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=memory_footprint, ylims=[0, 20], y_ticks = [0, 5, 10, 15, 20], x_label="Proportion Updated", y_label="Memory (MB)", grouping=proportion_updated, filename="just-edge-updates-memory")
