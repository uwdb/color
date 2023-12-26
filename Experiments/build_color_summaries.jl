function build_experiments(experiment_params_list::Vector{ExperimentParams})
    @sync @distributed for experiment_params in shuffle(experiment_params_list)
        build_times = [("Dataset", "Partitioner", "NumColors",  "BuildPhase", "BuildTime", "MemoryFootprint")]
        dataset = experiment_params.dataset
        summary_params = experiment_params.summary_params
        data = load_dataset(dataset)
        vertices_to_add = nv(data.graph) - round(convert(Float64, nv(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_updated))
        cloned_data = DataGraph(convert(Int64, vertices_to_add))
        vertices_for_later = []
        old_to_new_node_mapping = Dict()
        edges_for_later = []
        current_node = vertices_to_add
        # chooses a subset of the graph to be preloaded
        if (convert(Float64, experiment_params.summary_params.proportion_updated) > 0.0)
            graph_vertices = collect(vertices(data.graph))
            shuffle(graph_vertices)
            graph_edges = collect(edges(data.graph)) # these store the connections of the original nodes...
            for vertex in graph_vertices
                if (vertices_to_add > 0)
                    add_labeled_node!(cloned_data, cloned_data.vertex_labels[vertex])
                    old_to_new_node_mapping[vertex] = vertex
                    vertices_to_add -= 1
                else
                    push!(vertices_for_later, vertex)
                    old_to_new_node_mapping[vertex] = current_node
                    current_node += 1
                end
            end
            for edge in graph_edges
                if !(src(edge) in keys(vertices_for_later)) && !(dst(edge) in keys(vertices_for_later))
                    add_labeled_edge!(cloned_data, (src(edge), dst(edge)), only(data.edge_labels[(src(edge), dst(edge))]))
                else
                    push!(edges_for_later, edge)
                end
            end
        end
        summary_name = params_to_summary_filename(experiment_params)
        summary_file_location = "Experiments/SerializedSummaries/" * summary_name
        println("Building Color Summary: ", summary_name)
        timing_vec = Float64[]
        results = @timed generate_color_summary((experiment_params.summary_params.proportion_updated > 0) ? cloned_data : data, summary_params; verbose=1, timing_vec=timing_vec)
        current_summary = results.value
        # updates the remaining portion of the graph
        if (experiment_params.summary_params.proportion_updated > 0)
            for vertex in vertices_for_later
                add_summary_node!(current_summary, get(data.vertex_labels, vertex, []), vertex + 1) # need to double check if you do add one here
            end
            for edge in edges_for_later
                add_summary_edge!(current_summary, old_to_new_node_mapping[src(edge)], old_to_new_node_mapping[dst(edge)], get(data.edge_labels, (old_to_new_node_mapping[src(edge)], old_to_new_node_mapping[dst(edge)]), []))
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
