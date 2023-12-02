function build_experiments(experiment_params_list::Vector{ExperimentParams})
    for experiment_params in experiment_params_list
        build_times = [("Dataset", "Partitioner", "NumColors",  "BuildPhase", "BuildTime", "MemoryFootprint")]
        dataset = experiment_params.dataset
        summary_params = experiment_params.summary_params
        data = load_dataset(dataset)
        cloned_data = DataGraph(nv(data.graph))
        remaining_edges = []
        if (experiment_params.summary_params.proportion_not_updated < 1.0)
            cloned_data.vertex_labels = data.vertex_labels
            graph_edges = collect(edges(data.graph))
            # edges_to_add = (length(graph_edges) * experiment_params.summary_params.proportion_not_updated)
            for edge in graph_edges
                if (rand() < experiment_params.summary_params.proportion_not_updated)
                    add_labeled_edge!(cloned_data, (src(edge), dst(edge)), only(data.edge_labels[(src(edge), dst(edge))]))
                    # edges_to_add -= 1
                else
                    push!(remaining_edges, edge)
                end
            end
        end
        summary_name = params_to_summary_filename(experiment_params)
        summary_file_location = "Experiments/SerializedSummaries/" * summary_name
        println("Building Color Summary: ", summary_name)
        timing_vec = Float64[]
        results = @timed generate_color_summary((experiment_params.summary_params.proportion_not_updated < 1.0) ? cloned_data : data, summary_params; verbose=1, timing_vec=timing_vec)
        current_summary = results.value
        if (experiment_params.summary_params.proportion_not_updated < 1.0)
            for edge in remaining_edges
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
end
