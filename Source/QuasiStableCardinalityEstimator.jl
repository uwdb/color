# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
include("PropertyGraph.jl")
include("ExactSizeCalculator.jl")
include("ColorSummary.jl")
using AutoHashEquals
using Graphs
using QuasiStableColors
using Probably
using StatsBase
using Plots, GraphRecipes

# The following two functions sum over all paths which have the same color assigned to a particular node in the query graph.
# Equivalently, they perform a groupby on all other nodes of the query graph. The goal of this is to prevent
# an exponential growth in the number of paths through the lifted color graph. However, we can only remove query nodes whose
# edges have already been processed.
function sum_over_node!(partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}}, current_query_nodes, node_to_remove)
    nodeIdx = 1
    for node in current_query_nodes
        if node == node_to_remove
            break
        end
        nodeIdx += 1
    end
    new_partial_paths::Dict{Vector{Int}, Union{Vector{Float64}, Int}} = Dict()
    for path_and_bounds in partial_paths
        path = path_and_bounds[1]
        bounds = path_and_bounds[2]
        new_path = copy(path)
        deleteat!(new_path, nodeIdx)
        if !haskey(new_partial_paths, new_path)
            new_partial_paths[new_path] = bounds
        else
            new_partial_paths[new_path] = new_partial_paths[new_path] .+ bounds
        end
    end
    deleteat!(current_query_nodes, nodeIdx)
    empty!(partial_paths)
    for path in keys(new_partial_paths)
        push!(partial_paths, (path, new_partial_paths[path]))
    end
end

function sample_paths(partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}}, num_samples::Int)
    # partial_path[x] = (color path, bounds)
    partial_paths = [x for x  in partial_paths if x[2][2] > 0]

    # if we want to sample more paths than there are existing, then just return the original partial paths
    if (num_samples > length(partial_paths))
        return partial_paths
    end


    # sum up all of the bounds
    overall_bounds_sum::Vector{Float64} = [0,0,0]
    for path_and_bounds in partial_paths
        overall_bounds_sum = overall_bounds_sum .+ path_and_bounds[2]
    end

    # choose a sample of the paths
    sample_weights = [x[2][2] for x in partial_paths]
    sample_weights =AnalyticWeights(sample_weights./sum(sample_weights))
    path_samples::Vector{Tuple{Vector{Int}, Vector{Float64}}} = sample(partial_paths, sample_weights,  num_samples; replace=false)

    # sum up the sampled bounds
    sampled_bounds_sum::Vector{Float64} = [0,0,0]
    for path_and_bounds in path_samples
        sampled_bounds_sum = sampled_bounds_sum .+ path_and_bounds[2]
    end

    # get the difference between the overall and sampled bound sum_over_finished_query_nodes
    bound_diff::Vector{Float64} = overall_bounds_sum .- sampled_bounds_sum

    # for each sampled path...
    for i in eachindex(path_samples)
        # figure out what fraction of the sampled bounds is in the current Bounds
        # higher bounds will have more weight redistributed to them
        bound_fractions = path_samples[i][2] ./ sampled_bounds_sum
        # use that fraction of the difference (i.e. the removed path weights) and add it to the partial path
        redistributed_weights = bound_fractions .* bound_diff
        path_samples[i][2][1] = path_samples[i][2][1] + redistributed_weights[1]
        path_samples[i][2][2] = path_samples[i][2][2] + redistributed_weights[2]
        path_samples[i][2][3] = path_samples[i][2][3] + redistributed_weights[3]
        # path_samples[i][2] = path_samples[i][2] .+ redistributed_weights

        # not sure if I'm doing this right
        # if the weights go in terms of min/avg/max
        # if you combine two paths
        # what does the min/avg/max become?
        # does it become the min/avg/max of the two
    end

    return path_samples
end

function handle_extra_edges!(query::QueryGraph, summary::ColorSummary, partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}},
                                current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int,Int}}, usingStoredStats::Bool)
    # To account for cyclic queries, we check whether there are any remaining edges that have not
    # been processed. If so, we set the lower bound to 0, reduce the average estimate accordingly, and leave
    # the upper bound unchanged.
    remaining_edges = []
    for edge in edges(query.graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges) &&
                 (src(edge) in current_query_nodes && dst(edge) in current_query_nodes)
            push!(remaining_edges, (src(edge), dst(edge)))
            push!(visited_query_edges, (src(edge), dst(edge)))
        end
    end

    # scale down the average if there are remaining non-tree-edges
    for edge in remaining_edges
        parent_node_idx = indexin(edge[1], current_query_nodes)
        new_node_idx = indexin(edge[2], current_query_nodes)
        child_label = only(query.vertex_labels[edge[2]])
        edge_label = only(query.edge_labels[(edge[1],edge[2])])
        for i  in range(1, length(partial_paths))
            path = partial_paths[i][1]
            bounds = partial_paths[i][2]
            lower = 0
            average = only(bounds[2])
            parent_color = only(path[parent_node_idx])
            child_color = only(path[new_node_idx][1])
            current_colors::StartEndColorPair = [child_color, parent_color]
            # don't have to check data label because these nodes are already in the
            # partial path, so we have already ensured that the colors are appropriate
            probability_of_edge = 0
            if (haskey(summary.edge_avg_out_deg, edge_label)
                    && haskey(summary.edge_avg_out_deg[edge_label], child_label) # so we know that the child label is not appearing in the edge label table...
                        && haskey(summary.edge_avg_out_deg[edge_label][child_label], parent_color)
                            && haskey(summary.edge_avg_out_deg[edge_label][child_label][parent_color], child_color))
                if usingStoredStats
                    # we flip this because the matching graph finds the path between two nodes,
                    # where the last node is the start of the closing edge
                    path_graph = get_matching_graph(edge[2], edge[1], query)
                    path_bools = convert_path_graph_to_bools(path_graph)
                    current_cycle_description = CyclePathAndColors(path_bools, current_colors)
                    default_cycle_description = CyclePathAndColors(path_bools, (-1, -1))
                    if haskey(summary.cycle_probabilities, current_cycle_description)
                        probability_of_edge = summary.cycle_probabilities[current_cycle_description]
                    elseif haskey(summary.cycle_probabilities, default_cycle_description)
                        probability_of_edge = summary.cycle_probabilities[default_cycle_description]
                    end
                else
                    probability_of_edge = get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary)
                end
            end
            average *= probability_of_edge
            upper = only(bounds[3])
            partial_paths[i] = (path, [lower, average, upper])
        end
    end
end

function sum_over_finished_query_nodes!(query::QueryGraph, partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}},
                                            current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int, Int}})
    prev_query_nodes = copy(current_query_nodes)
    for node in prev_query_nodes
        has_living_edges = false
        for neighbor in outneighbors(query.graph, node)
            if !((node, neighbor) in visited_query_edges)
                has_living_edges = true
            end
        end
        for neighbor in inneighbors(query.graph, node)
            if !((neighbor, node) in visited_query_edges)
                has_living_edges = true
            end
        end
        if ! has_living_edges
            sum_over_node!(partial_paths, current_query_nodes, node)
        end
    end
end

function get_min_width_node_order(g::DiGraph)
    if (nv(g) < 10)
        nodes_processed = 1
        partial_orders::Dict{Set{Int64}, Tuple{Vector{Int64}, Int64}} = Dict(Set([x]) =>([x], 0) for x in vertices(g))
        while nodes_processed < nv(g)
            for node_set in keys(partial_orders)
                best_order, best_width = partial_orders[node_set]
                neighbor_nodes = Set()
                for existing_node in node_set
                    for neighbor in all_neighbors(g, existing_node)
                        if !(neighbor in node_set)
                            push!(neighbor_nodes, neighbor)
                        end
                    end
                end
                for next_node in neighbor_nodes
                    new_node_set::Set{Int32} = Set([next_node])
                    new_node_set = union(new_node_set, node_set)
                    new_order = []
                    copy!(new_order, best_order)
                    push!(new_order, next_node)
                    new_width = 0
                    for v in new_order
                        if ! all([x in new_order for x in all_neighbors(g, v)])
                            new_width += 1
                        end
                    end
                    new_width = max(best_width, new_width)
                    if haskey(partial_orders, new_node_set)
                        if partial_orders[new_node_set][2] > new_width
                            partial_orders[new_node_set] = (new_order, new_width)
                        end
                    else
                        partial_orders[new_node_set] = (new_order, new_width)
                    end
                end
            end
            nodes_processed += 1
            for node_set in keys(partial_orders)
                if length(node_set) < nodes_processed
                    delete!(partial_orders, node_set)
                end
            end
        end
        min_width = minimum([x[2] for x in values(partial_orders)])
        for node_order_and_width in values(partial_orders)
            if node_order_and_width[2] == min_width
                return node_order_and_width[1]
            end
        end
    else
        min_width = nv(g)
        min_order = []
        for starting_node in vertices(g)
            max_width = 0
            visited_nodes = [starting_node]
            while length(visited_nodes) < nv(g)
                new_width = nv(g)
                next_node = -1
                for potential_node in vertices(g)
                    if potential_node in visited_nodes || !any([x in all_neighbors(g, potential_node) for x in visited_nodes])
                        continue
                    end
                    potential_visited_nodes = []
                    copy!(potential_visited_nodes, visited_nodes)
                    push!(potential_visited_nodes, potential_node)
                    potential_num_active_nodes = 0
                    for v in potential_visited_nodes
                        if ! all([x in potential_visited_nodes for x in all_neighbors(g, v)])
                            potential_num_active_nodes += 1
                        end
                    end
                    if potential_num_active_nodes <= new_width
                        next_node = potential_node
                        new_width = potential_num_active_nodes
                    end
                end
                push!(visited_nodes, next_node)
                max_width = max(max_width, new_width)
            end
            if max_width <= min_width
                min_order = visited_nodes
                min_width = max_width
            end
        end
        return min_order
    end
end

function get_cardinality_bounds(query::QueryGraph, summary::ColorSummary; max_partial_paths = nothing, use_partial_sums = true, verbose = false, usingStoredStats = false, include_cycles = true)
    node_order = get_min_width_node_order(query.graph) #spanning tree to cut out cycles
    if verbose
        println("Node Order:", node_order)
    end
    # Because the label is implied by the color -> query_graph_vertex mapping stored in current_query_nodes,
    # we don't have to keep the label in the partial paths object.
    partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}} = [] # each tuple contains a pairing of color paths -> bounds
    visited_query_edges::Vector{Tuple{Int,Int}} = []
    current_query_nodes::Vector{Int} = []

    old_node = popfirst!(node_order)
    parent_label = query.vertex_labels[old_node][1]
    parent_data_labels = get_data_label(query, old_node)
    push!(current_query_nodes, old_node)
    # Initialize partial_paths with all possible starting color/vertex possibilities.
    for color in keys(summary.color_label_cardinality)
        # Only use the parent label.
        if (haskey(summary.color_label_cardinality[color], parent_label))
            # if the parent has a specified data label, only use colors that the filters approve
            data_label_is_in_color = false
            for data_label in parent_data_labels
                if data_label == -1
                    data_label_is_in_color = true
                    continue
                end
                if data_label in summary.color_filters[color]
                    data_label_is_in_color = true
                end
            end
            if (data_label_is_in_color)
                push!(partial_paths, ([color], [summary.color_label_cardinality[color][parent_label],
                                                            summary.color_label_cardinality[color][parent_label],
                                                            summary.color_label_cardinality[color][parent_label]]))
            end
        end
    end

    new_node = old_node
    while length(node_order) > 0
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        if use_partial_sums
            sum_over_finished_query_nodes!(query, partial_paths, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end

        # Get the next child from the node order.
        new_node = popfirst!(node_order)
        parent_idx = 0
        outEdge = false
        for neighbor in all_neighbors(query.graph, new_node)
            if neighbor in current_query_nodes
                old_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
                if old_node in inneighbors(query.graph, new_node)
                    outEdge = true
                end
                break
            end
        end
        # Push the current edge and nodes to the visited lists.
        if outEdge
            push!(visited_query_edges, (old_node, new_node))
        else
            push!(visited_query_edges, (new_node, old_node))
        end
        push!(current_query_nodes, new_node)

        # Get the appropriate labels, the query only uses one label per vertex/node.
        edge_label = 0
        if outEdge
            edge_label = only(query.edge_labels[(old_node,new_node)])
        else
            edge_label =  only(query.edge_labels[(new_node,old_node)])
        end
        new_label = only(query.vertex_labels[new_node])
        new_data_labels = get_data_label(query, new_node)

        # Update the partial paths using the parent-child combo that comes next from the query.
        new_partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}} = []
        if (max_partial_paths !== nothing)
            if (length(partial_paths) > max_partial_paths)
                partial_paths = sample_paths(partial_paths, max_partial_paths)
            end
        end
        for path_and_bounds in partial_paths
            path = path_and_bounds[1]
            running_bounds = path_and_bounds[2]
            old_color = only(path[parent_idx])
            # Account for colors with no outgoing children.
            if outEdge && haskey(summary.edge_avg_out_deg, edge_label) &&
                haskey(summary.edge_avg_out_deg[edge_label], new_label) &&
                haskey(summary.edge_avg_out_deg[edge_label][new_label], old_color)
                for new_color in keys(summary.edge_avg_out_deg[edge_label][new_label][old_color])
                    # revamp the logic to use a set of labels rather than just one

                    # check if the data label(s) are in the color
                    data_label_in_color = false
                    for data_label in new_data_labels
                        if data_label == -1
                            data_label_in_color = true
                            continue
                        end
                        if data_label in summary.color_filters[new_color]
                            data_label_in_color = true
                        end
                    end
                    if !data_label_in_color
                        continue
                    end

                    new_path = copy(path)
                    push!(new_path, new_color)
                    new_bounds = [running_bounds[1]*summary.edge_min_out_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[2]*summary.edge_avg_out_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[3]*summary.edge_max_out_deg[edge_label][new_label][old_color][new_color],
                                    ]
                    if (new_data_labels != [-1])
                        # we have already confirmed that the data label is in the color, but if the data label isn't -1
                        # then we need to scale down the result since we only want to consider one of the many nodes in the new color
                        new_bounds[2] = new_bounds[2] / summary.color_label_cardinality[new_color][new_label]
                        # we also need to set the minimum to 0 but keep the maximum the same
                        new_bounds[1] = 0
                    end
                    push!(new_partial_paths, (new_path, new_bounds))
                end
            elseif !outEdge && haskey(summary.edge_avg_in_deg, edge_label) &&
                    haskey(summary.edge_avg_in_deg[edge_label], new_label) &&
                    haskey(summary.edge_avg_in_deg[edge_label][new_label], old_color)
                for new_color in keys(summary.edge_avg_in_deg[edge_label][new_label][old_color])
                    data_label_in_color = false
                    for data_label in new_data_labels
                        if data_label == -1
                            data_label_in_color = true
                            continue
                        end
                        if data_label in summary.color_filters[new_color]
                            data_label_in_color = true
                        end
                    end
                    if !data_label_in_color
                        continue
                    end
                    new_path = copy(path)
                    push!(new_path, new_color)
                    new_bounds = [running_bounds[1]*summary.edge_min_in_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[2]*summary.edge_avg_in_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[3]*summary.edge_max_in_deg[edge_label][new_label][old_color][new_color],
                                    ]
                    if (new_data_labels != [-1])
                        # we have already confirmed that the data label is in the color, but if the data label isn't -1
                        # then we need to scale down the result since we only want to consider one of the many nodes in the new color
                        new_bounds[2] = new_bounds[2] / summary.color_label_cardinality[new_color][new_label]
                        # we also need to set the minimum to 0 but keep the maximum the same
                        new_bounds[1] = 0
                    end
                    push!(new_partial_paths, (new_path, new_bounds))
                end
            end
        end
        partial_paths = new_partial_paths
        if (include_cycles)
            handle_extra_edges!(query, summary, partial_paths, current_query_nodes, visited_query_edges, usingStoredStats)
        end
    end

    # Sum over the calculated partial paths to get the final bounds.
    final_bounds = [0,0,0]
    for path_and_bounds in partial_paths
        final_bounds = final_bounds .+ path_and_bounds[2]
    end
    return final_bounds
end
