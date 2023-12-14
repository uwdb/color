using Plots.PlotMeasures
using Graphs
using Random
include("../Experiments.jl")

datasets::Vector{DATASET} = [aids]
# datasets::Vector{DATASET} = [aids, human, yeast, wordnet, youtube, dblp, patents]
# datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
max_cycles = 6
proportions_deleted = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
# To test deletion, we will add a random node / edge and then delete them...
# proportion_not_updated = 0.5

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, max_cycle_size=current_cycle, proportion_deleted=current_proportion) 
                                                    for current_dataset in datasets for current_cycle in 2:max_cycles for current_proportion in proportions_deleted]
println("started building")
    for experiment_params in experiment_params_list
        build_times = [("Dataset", "Partitioner", "NumColors",  "BuildPhase", "BuildTime", "MemoryFootprint")]
        dataset = experiment_params.dataset
        summary_params = experiment_params.summary_params
        data = load_dataset(dataset)
        vertices_to_add = round(convert(Float64, nv(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_not_updated))
        cloned_data = DataGraph(convert(Int64, vertices_to_add))
        vertices_for_later = []
        old_to_new_node_mapping = Dict()
        edges_for_later = []
        current_node = vertices_to_add
        if (convert(Float64, experiment_params.summary_params.proportion_not_updated) < 1.0)
            graph_vertices = collect(vertices(data.graph))
            shuffle(graph_vertices)
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
        results = @timed generate_color_summary((experiment_params.summary_params.proportion_not_updated < 1.0) ? cloned_data : data, summary_params; verbose=1, timing_vec=timing_vec)
        current_summary = results.value
        if (experiment_params.summary_params.proportion_not_updated < 1.0)
            for vertex in vertices_for_later
                add_summary_node!(current_summary, get(data.vertex_labels, vertex, []), vertex + 1) # need to double check if you do add one here
            end
            for edge in edges_for_later
                add_summary_edge!(current_summary, old_to_new_node_mapping[src(edge)], old_to_new_node_mapping[dst(edge)], get(data.edge_labels, (old_to_new_node_mapping[src(edge)], old_to_new_node_mapping[dst(edge)]), []))
            end
        end
        # handle deletion testing here:
        temp_vertices = Dict() # map temp nodes to their labels
        temp_edges = Dict() # map temp edges to their labels
        # first, keep track of which nodes you want to delete and then delete them
        temp_vertices_to_delete = round(convert(Float64, nv(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_deleted))
        # for each node, delete all edges related to them (and keep track of it), then remove the node
        for i in 1:temp_vertices_to_delete
            current_node = rand(1:nv(cloned_data.graph))
            temp_vertices[current_node] = cloned_data.vertex_labels[current_node]
            for edge in edges(cloned_data.graph)
                if (src(edge) == current_node || dst(edge) == current_node)
                    if !(edge in keys(temp_edges))
                        temp_edges[edge] = cloned_data.edge_labels[(src(edge), dst(edge))]
                        remove_summary_edge!(current_summary, src(edge), dst(edge), [])
                    end
                end
            end
            delete_summary_node!(current_summary, temp_vertices[current_node], current_node)
        end
        # now, go through all of the removed vertices and add them, then add all of the removed edges
        for vertex in keys(temp_vertices)
            add_summary_node!(current_summary, temp_vertices[vertex], vertex)
        end
        for edge in keys(temp_edges)
            add_summary_edge!(current_summary, src(edge), dst(edge), temp_edges[edge])
        end
        # if successful, the results should be comparable to the regular results for the queries
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
# compare how overall accuracy is affected by summary updates
# graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=proportion_not_updated, filename="overall-accuracy-and-updates")
# compare how cycle stat accuracies are affected by summary updates
# graph_grouped_box_plot(experiment_params_list, x_type=proportion_deleted, y_type=estimate_error, x_label="proportion added then deleted", y_label="accuracy", grouping=cycle_size, filename="deletion-experiment")
graph_grouped_bar_plot(experiment_params_list, x_type=proportion_deleted, y_type=build_time, x_label="Proportion Deleted Then Updated", y_label="Build Time", grouping=dataset, filename="delete-build-aids")
graph_grouped_box_plot(experiment_params_list, x_type=proportion_deleted, y_type=estimate_error, x_label="Proportion Deleted Then Updated", y_label="Estimate Error", grouping=dataset, filename="delete-error-aids")
graph_grouped_bar_plot(experiment_params_list, x_type=proportion_deleted, y_type=runtime, y_lims=[0, 0.02], x_label="Proportion Deleted Then Updated", y_label="Runtime", grouping=dataset, filename="delete-runtime-aids")
graph_grouped_bar_plot(experiment_params_list, x_type=proportion_deleted, y_type=memory_footprint, y_lims=[0, 6], x_label="Proportion Deleted Then Updated", y_label="Memory Footprint", grouping=dataset, filename="delete-memory-aids")