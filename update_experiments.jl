using Plots.PlotMeasures
using Graphs
include("Experiments/Experiments.jl")

datasets::Vector{DATASET} = [human]
# datasets::Vector{DATASET} = [aids, human, yeast, wordnet, youtube, dblp, patents]
# datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
max_cycles = 6
proportions_not_updated = [0, 0.2, 0.4, 0.6, 0.8, 1]
# proportions_not_updated = [1.0, 0.8]


for proportion_not_updated in proportions_not_updated
# only do this with the cloned stuff... we can repeat the same experiments without cloning in the normal version of the code
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_size) for current_dataset in datasets for current_size in 2:max_cycles]
println("started building")
for experiment_params in experiment_params_list
    build_times = [("Dataset", "Partitioner", "NumColors", "BuildTime", "MemoryFootprint")]
    dataset = experiment_params.dataset
    summary_params = experiment_params.summary_params
    data = load_dataset(dataset)
    cloned_data = DataGraph(nv(data.graph)) # remember to add back in the vertex labels
    cloned_data.vertex_labels = data.vertex_labels
    graph_edges = collect(edges(data.graph))
    println("length: ", length(graph_edges))
    edges_to_add = (length(graph_edges) * proportion_not_updated)
    println("edges to add: ", length(graph_edges) * proportion_not_updated)
    remaining_edges = []
    for edge in graph_edges
        if (edges_to_add > 0)
            add_labeled_edge!(cloned_data, (src(edge), dst(edge)), only(data.edge_labels[(src(edge), dst(edge))]))
            edges_to_add -= 1
            # update_edge_labels!(cloned_data, (src(edge), dst(edge)), data.edge_labels[(src(edge), dst(edge))])
        else
            push!(remaining_edges, edge)
        end
    end
    summary_name = (params_to_summary_filename(experiment_params) * "with" * string(proportion_not_updated) * "updates")
    summary_file_location = "Experiments/SerializedSummaries/" * summary_name
    println("Building Color Summary: ", summary_name)
    results= @timed generate_color_summary(cloned_data, summary_params; verbose=1, detailed_cycles=false)    # println("normal time: ", normal_results.time)
    current_summary = results.value
    for edge in remaining_edges
        # println("ADDING SUMMARY EDGES")
        add_summary_edge!(current_summary, src(edge), dst(edge), get(data.edge_labels, (src(edge), dst(edge)), []))
    end
    summary_size = Base.summarysize(current_summary)
    serialize(summary_file_location, current_summary)
    push!(build_times, (string(dataset),
                         string(summary_params.partitioner),
                         string(summary_params.num_colors),
                         string(results.time),
                         string(summary_size)))
    results_filename = params_to_results_filename(experiment_params)
    result_file_location = "Experiments/Results/Build_" * results_filename
    writedlm(result_file_location, build_times, ",")
end
println("started estimating")
for experiment_params in experiment_params_list
    dataset = experiment_params.dataset
    all_queries = load_querysets([dataset]; require_true_cardinality = true)
    summary_name = (params_to_summary_filename(experiment_params) * "with" * string(proportion_not_updated) * "updates")
    summary_file_location = "Experiments/SerializedSummaries/" * summary_name
    !isfile(summary_file_location) && error("The summary has not been built yet! \n Attempted File Location: $(summary_file_location)")
    summary::ColorSummary = deserialize(summary_file_location)
    experiment_results = []
    push!(experiment_results, ("UpperBound", "Estimate", "LowerBound", "TrueCard", "EstimationTime", "QueryType"))
    for i in 1:length(all_queries[dataset])
        query = all_queries[dataset][i].query
        query_path = all_queries[dataset][i].query_path
        exact_size = all_queries[dataset][i].exact_size
        results = @timed get_cardinality_bounds(query, summary;
                            max_partial_paths = experiment_params.inference_max_paths,
                            use_partial_sums=experiment_params.use_partial_sums, usingStoredStats=true,
                            sampling_strategy=experiment_params.sampling_strategy,
                            only_shortest_path_cycle= experiment_params.only_shortest_path_cycle)
        upper_bound = results.value[3]
        estimate = max(1, results.value[2])
        lower_bound = results.value[1]
        estimate_time = results.time
        query_type = all_queries[dataset][i].query_type
        push!(experiment_results, (upper_bound, estimate, lower_bound, exact_size, estimate_time, query_type))
    end
    results_file_location = "Experiments/Results/Estimation_"  * params_to_results_filename(experiment_params) * "with" * string(proportion_not_updated) * "updates"
    writedlm(results_file_location, experiment_results, ",")
end
println("started graphing")
# graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=cycle_size, filename="detailed-cycle-experiment")
x_type::GROUP=dataset
y_type::VALUE=estimate_error
grouping::GROUP=cycle_size    
x_values = []
    y_values = []
    groups = []
    for experiment_params in experiment_params_list
        # load the results
        results_filename = params_to_results_filename(experiment_params)
        results_path = "Experiments/Results/Estimation_" * results_filename * "with" * string(proportion_not_updated) * "updates"
        # println("results path: ", results_path)
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)

        # get the x_value and grouping (same for all results in this experiment param)

        # keep track of the data points
        for i in 1:nrow(results_df)
            current_x = x_type == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, x_type)
            current_group = grouping == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, grouping)
            current_y = 0
            if y_type == estimate_error
                current_y = results_df[i, :Estimate] / results_df[i, :TrueCard]
            else # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, current_x)
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = groupedboxplot(x_values, y_values, group = groups, yscale =:log10,
                            ylims=[10^-13, 10^11], yticks=[10^-10, 10^-5, 10^-2, 1, 10^2, 10^5, 10^10],
                            legend = :outertopleft, size = (1000, 600))
    xlabel!(gbplot, "Dataset")
    ylabel!(gbplot, "Accuracy")
    plotname = "cycleswith" * string(proportion_not_updated) * "withcycleupdateshuman" * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)

# first, take the data graph and collect its list of edges
# then, create a new data graph with the same number of nodes
# from the original graph, take the list of edges
# for a portion of the edges, put them into the new graph directly
# for the remaining edges, add them as "summary edges" using the summary
# then, run a query, using the "real" graph and the "updated" graph
end