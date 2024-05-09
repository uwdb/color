using Plots.PlotMeasures
using Graphs
using Random
include("../Experiments.jl")

# The goal of this file is to determine the performance of the node and edge update logic.
# We create a summary using only some of the edges (and their involved nodes) in the graph and then update the summary
# with the remaining edges (and any unadded nodes), repeating for different proportions of the graph to update with.
# Then, we can check metrics like relative error and see how it performs with more updates.

datasets::Vector{DATASET} = [aids]
proportions_updated = [0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, proportion_updated=current_proportion)
                                                    for current_dataset in datasets for current_proportion in proportions_updated]
println("started building")
shuffled_edges = Dict()
for experiment_params in experiment_params_list
    build_times = [("Dataset", "Partitioner", "NumColors", "BuildPhase", "BuildTime", "MemoryFootprint", "VertexUpdateTime", "EdgeUpdateTime", "NumVertices", "NumEdges")]
    dataset = experiment_params.dataset
    summary_params = experiment_params.summary_params
    println("Current proportion: ", summary_params.proportion_updated)

    data = load_dataset(dataset)
    println("Finished loading dataset")
    cloned_data = DataGraph(0)
    edges_for_later = []
    added_nodes = Set()
    nodes_mapping = Dict() # map the current node in the cloned graph to its node id in the original graph
    # whenever we add edges, we have to refer back to the node mapping
    shuffled_edges[dataset] = shuffle!(collect(edges(data.graph)))
    num_start_edges = ne(data.graph) - round(convert(Float64, ne(data.graph)) * convert(Float64, experiment_params.summary_params.proportion_updated))
    num_start_edges = convert(Int, num_start_edges)
    # grab the list of edges
    # create a graph with the number of unique edges making of the start/end of the edges
    # make a mapping of cloned graph node id => original node id
    # update the vertex labels according to the original graph's labels for those nodes
    # for each edge, go through and add them to the graph
    edges_to_add = shuffled_edges[dataset][1:num_start_edges]
    nodes_to_add = Set()
    for edge in edges_to_add
        push!(nodes_to_add, src(edge))
        push!(nodes_to_add, dst(edge))
    end
    cloned_data = DataGraph(length(nodes_to_add))
    vertex_mapping = Dict()
    vertices_for_later = []
    println("Creating initial graph")
    # vertices start at 1, data labels start at 0...
    vertex_in_clone = 0
    for edge in shuffled_edges[dataset]
        current_vertices = (src(edge), dst(edge))
        if edge in edges_to_add
            # for each node in the edge, track its mapping and update its label
            for current_vertex in current_vertices
                if !(current_vertex in added_nodes)
                    vertex_in_clone += 1
                    vertex_mapping[current_vertex] = vertex_in_clone
                    update_node_labels!(cloned_data, vertex_mapping[current_vertex], data.vertex_labels[current_vertex])
                end
                push!(added_nodes, current_vertex)
            end
            # now explicitly add the edge to the graph
            mapped_start = vertex_mapping[current_vertices[1]]
            mapped_end = vertex_mapping[current_vertices[2]]
            for label in data.edge_labels[current_vertices]
                add_labeled_edge!(cloned_data, (mapped_start, mapped_end), label)
            end
        else
            # save the edge (and its vertices) to be processed later
            push!(edges_for_later, edge)
        end
    end
    summary_name = params_to_summary_filename(experiment_params)
    summary_file_location = "Experiments/SerializedSummaries/" * summary_name
    println("Building Color Summary: ", summary_name)
    timing_vec = Float64[]
    results = @timed generate_color_summary((experiment_params.summary_params.proportion_updated > 0) ? cloned_data : data, summary_params; verbose=1, timing_vec=timing_vec)
    current_summary = results.value
    vertex_update_time = 0
    edge_update_time = 0
    num_edges = length(edges_for_later)
    num_vertices = nv(data.graph) - length(added_nodes)
    println("Updating graph as necessary")
    if (experiment_params.summary_params.proportion_updated > 0)
        for edge in edges_for_later
            # first add the summary nodes to the graph if need
            current_vertices = (src(edge), dst(edge))
            for current_vertex in current_vertices
                # update the mapping if the nodes haven't already been added
                if !(current_vertex in added_nodes)
                    vertex_in_clone += 1
                    vertex_mapping[current_vertex] = vertex_in_clone
                    # now add a summary node_label
                    vertex_update_time += @elapsed add_summary_node!(current_summary, data.vertex_labels[current_vertex], vertex_in_clone)
                end
                push!(added_nodes, current_vertex)
            end
            # now add the summary edge to the graph
            mapped_start = vertex_mapping[current_vertices[1]]
            mapped_end = vertex_mapping[current_vertices[2]]
            edge_update_time += @elapsed add_summary_edge!(current_summary, mapped_start, mapped_end, data.edge_labels[current_vertices])
        end
    end
    
    println("Serializing")
    summary_size = Base.summarysize(current_summary)
    serialize(summary_file_location, current_summary)
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "FullTime",
        string(results.time),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "Coloring",
        string(timing_vec[1]),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "CycleCounting",
        string(timing_vec[2]),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "BloomFilter",
        string(timing_vec[3]),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "CardinalityCounting",
        string(timing_vec[4]),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    push!(build_times, (string(dataset),
        string(summary_params.partitioning_scheme),
        string(summary_params.num_colors),
        "EdgeStats",
        string(timing_vec[5]),
        string(summary_size),
        string(vertex_update_time),
        string(edge_update_time),
        string(num_vertices),
        string(num_edges)))
    results_filename = params_to_results_filename(experiment_params)
    result_file_location = "Experiments/Results/Build_" * results_filename
    writedlm(result_file_location, build_times, ",")
end

println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_box_plot(experiment_params_list, dimensions=(600,400), x_type=proportion_updated, y_type=estimate_error, ylims=[10e-10, 10e20], y_ticks=[10^-5, 1, 10^5, 10^10, 10^15], x_label="Proportion Updated (AIDS)", y_label="Relative Error (log\$_{10}\$)", grouping=proportion_updated, legend_pos=:outerright, filename="fig_11")
