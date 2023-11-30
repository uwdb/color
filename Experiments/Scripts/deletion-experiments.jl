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

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_cycle, proportion_deleted=current_proportion) 
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
        temp_vertices = []
        temp_edges = []
        # first, add and keep track of a specific number of nodes as you add them to the graph
        temp_vertices_to_add = round(convert(Float64, nv(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_deleted))
        # for each node, choose an already-existing node and add and keep track of an edge between them.
        for i in 1:temp_vertices_to_add
            current_node = nv(cloned_data.graph) + i
            add_summary_node!(current_summary, [], current_node)
            other_node = rand(0:nv(cloned_data.graph))
            add_summary_edge!(current_summary, current_node, other_node, [])
            push!(temp_vertices, current_node)
            push!(temp_edges, (current_node, other_node))
        end
        # now, go through all of the added edges and delete them, then go through all of the added nodes and delete them
        for edge in temp_edges
            remove_summary_edge!(current_summary, edge[1], edge[2], [])
        end
        for vertex in temp_vertices
            delete_summary_node!(current_summary, [], current_node)
        end
        # if successful, the results should be comparable to the regular results for the queries
        summary_size = Base.summarysize(current_summary)
        serialize(summary_file_location, current_summary)
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "FullTime",
                             string(results.time),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "Coloring",
                             string(timing_vec[1]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "CycleCounting",
                             string(timing_vec[2]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "BloomFilter",
                             string(timing_vec[3]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "CardinalityCounting",
                             string(timing_vec[4]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
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
graph_grouped_box_plot(experiment_params_list, x_type=proportion_deleted, y_type=estimate_error, x_label="proportion added then deleted", y_label="accuracy", grouping=cycle_size, filename="deletion-experiment")