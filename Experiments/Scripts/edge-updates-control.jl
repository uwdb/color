using Plots.PlotMeasures
using Graphs
using Random
include("../Experiments.jl")

# The goal of this file is to verify that the edge update code is making meaningful adjustments
# to the stored edge statistics. To do this, we only partially build the graph and do not update
# with the remaining edges. Plots where the summaries were actually updated should demonstrate
# improved results.

datasets::Vector{DATASET} = [aids, human, wordnet, dblp]
proportions_updated = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

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
    # remove this part so that the edges aren't updated
    # if (experiment_params.summary_params.proportion_updated > 0)
    #     for edge in edges_for_later
    #         add_summary_edge!(current_summary, src(edge), dst(edge), get(data.edge_labels, (src(edge), dst(edge)), []))
    #     end
    # end
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
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=build_time, ylims=[0, 30], x_label="Proportion Updated", y_label="Build Time (s)", grouping=proportion_updated, filename="control-edge-updates-build")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, x_label="Proportion Updated", y_label="Estimate Error", grouping=proportion_updated, filename="control-edge-updates-error")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=runtime, ylims=[0, 0.6], x_label="Proportion Updated", y_label="Runtime (s)", grouping=proportion_updated, filename="control-edge-updates-runtime")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=memory_footprint, ylims=[0, 20], x_label="Proportion Updated", y_label="Memory Footprint (B)", grouping=proportion_updated, filename="control-edge-updates-memory")