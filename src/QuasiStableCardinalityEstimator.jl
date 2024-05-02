# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.

"""
Sums over all paths which have the same color assigned to a particular node in the query graph.
Equivalently, performs a groupby on all other nodes of the query graph. The goal of this is to prevent
an exponential growth in the number of paths through the lifted color graph. However, we can only remove query nodes whose
edges have already been processed. Returns an updated partial_paths and partial_weights.
# Arguments
- partial_paths::Matrix{Color} - the matrix of partial paths describing the current traversal over the lifted color graph.
- partial_weights::Vector{W} - the list of weights for each partial path.
- current_query_nodes - the list of query nodes, describing the order we process the query graph.
- node_to_remove - the specific node to sum over from the query.
"""
@inline function sum_over_node(partial_paths::Matrix{Color}, partial_weights::Vector{W}, current_query_nodes, node_to_remove) where W
    nodeIdx = 1
    for node in current_query_nodes
        if node == node_to_remove
            break
        end
        nodeIdx += 1
    end
    new_partial_paths::Dict{Vector{Color}, W} = Dict()
    for i in 1:size(partial_paths)[2]
        new_path = copy(partial_paths[:, i])
        deleteat!(new_path, nodeIdx)
        if !haskey(new_partial_paths, new_path)
            new_partial_paths[new_path] = partial_weights[i]
        else
            new_partial_paths[new_path] = sum_colorings(new_partial_paths[new_path], partial_weights[i])
        end
    end
    deleteat!(current_query_nodes, nodeIdx)
    partial_paths = zeros(Color, length(current_query_nodes), length(keys(new_partial_paths)))
    partial_weights = fill(W(0), length(keys(new_partial_paths)))

    path_idx = 1
    for path in keys(new_partial_paths)
        for i in eachindex(path)
            partial_paths[i, path_idx] = path[i]
        end
        partial_weights[path_idx] = new_partial_paths[path]
        path_idx += 1
    end
    return partial_paths, partial_weights
end

@enum SAMPLING_STRATEGY uniform weighted redistributive

"""
Samples the partial paths, returning a condensed list of partial paths. The partial weights are redistributed according
to the magnitude of the original weight of the output paths.
# Arguments
- partial_paths::Matrix{Color} - the matrix of partial paths describing the current traversal over the lifted color graph.
- partial_weights::Vector{W} - the list of weights for each partial path.
- num_samples::Int - the number of paths to sample, determining the size of the output.
- sampling_strategy::SAMPLING_STRATEGY - how to select samples from the partial paths. Will only check if this is set to "uniform" to sample uniformly.
                                         Otherwise, just selects samples by prioritizing those with higher weights.
"""
function sample_paths(partial_paths::Matrix{Color}, partial_weights::Vector{W}, num_samples::Int, sampling_strategy::SAMPLING_STRATEGY) where W
    # if we want to sample more paths than there are existing nonzero paths,
    # then just return the original partial paths
    new_partial_paths = zeros(Color, size(partial_paths))
    new_partial_weights = fill(W(0), size(partial_weights))
    new_path_idx = 1
    for i in eachindex(partial_weights)
        if get_count(partial_weights[i]) > 0
            new_partial_paths[:, new_path_idx] = partial_paths[:, i]
            new_partial_weights[new_path_idx] = partial_weights[i]
            new_path_idx += 1
        end
    end
    new_path_idx -= 1
    new_partial_paths = new_partial_paths[:, 1:new_path_idx]
    new_partial_weights = new_partial_weights[1:new_path_idx]

    if length(new_partial_weights) < num_samples
        return new_partial_paths, new_partial_weights
    end

    # sum up all of the bounds
    overall_bounds_sum = 0.0
    for w in new_partial_weights
        overall_bounds_sum += get_count(w)
    end

    # choose a sample of the paths
    sample_weights = [get_count(w) for w in new_partial_weights]
    sample_weights = AnalyticWeights(sample_weights ./ overall_bounds_sum)
    if sampling_strategy == uniform
        sample_weights = AnalyticWeights([1.0 for i in eachindex(new_partial_weights)] ./ length(new_partial_weights))
    end
    sample_indices::Vector{Int} = sample(1:length(new_partial_weights), sample_weights,  num_samples; replace=false)

    # sum up the sampled bounds
    sampled_bounds_sum = 0.0
    for idx in sample_indices
        sampled_bounds_sum += get_count(new_partial_weights[idx])
    end
    sampled_partial_paths = zeros(Color, size(new_partial_paths)[1], length(sample_indices))
    sampled_partial_weights = fill(W(0), length(sample_indices))

    for i in eachindex(sample_indices)
        idx = sample_indices[i]
        sampled_partial_paths[:, i] .= new_partial_paths[:, idx]
        inverse_sampling_probability = if sampling_strategy == redistributive
            # scale the weights so that their sum equals the input weight's sum
            overall_bounds_sum / sampled_bounds_sum
        else
            1.0 / (sample_weights[i] * num_samples)
        end
        sampled_partial_weights[i] = scale_coloring(new_partial_weights[idx], inverse_sampling_probability)
    end
    return sampled_partial_paths, sampled_partial_weights
end

"""
Creates a vector of paths from the start to destination node using a recursive depth-first search.
# Arguments
- visited::Set{Int} - a set of all the nodes that have already been processed during the path search.
- cur::Int - the starting node of the path.
- finish::Int - the destination node of the path.
- max_length::Int - the maximum length of paths to include in the search result.
- graph::SimpleGraph - the undirected version of the graph, to make path searches simpler.
- current_path::Vector{Int} - the current path we are traversing.
- simple_paths::Vector{Vector{Int}} - the list of resulting paths from the search.
"""
function get_simple_paths_dfs!(visited::Set{Int}, cur::Int, finish::Int, max_length::Int,
                                graph::SimpleGraph, current_path::Vector{Int},
                                simple_paths::Vector{Vector{Int}})
    length(current_path) > max_length && return
    cur in visited && return
    push!(visited, cur)
    push!(current_path, cur)
    if cur == finish
        push!(simple_paths, deepcopy(current_path))
        delete!(visited, cur)
        pop!(current_path)
        return
    end

    for next in all_neighbors(graph, cur)
        get_simple_paths_dfs!(visited, next, finish, max_length, graph, current_path,
                                simple_paths)
    end
    if length(current_path) > 0
        pop!(current_path)
    end
    delete!(visited, cur)
end

"""
Finds all the unique paths from the start to the end node in the given Query Graph.
Returns a list of paths, where each path is a Bool Vector representing the direction of the edge (true = forward, false = backward).
# Arguments
- start::Int - the starting node of the path in the Query Graph.
- finish::Int - the destination node of the path in the Query Graph.
- max_length::Int - the maximum length of the paths to return.
- query_graph::DiGraph - the Query Graph that the path belongs to.
- visited_edges::Vector{Tuple{Int,Int}} - a list of edges that have already been processed.
"""
function get_all_simple_path_bools(start::Int, finish::Int, max_length::Int,
                                    query_graph::DiGraph, visited_edges::Vector{Tuple{Int,Int}})
    # convert the graph to be undirected and only include the edges that have already been processed
    graph_copy = Graph(nv(query_graph))
    for edge in visited_edges
        add_edge!(graph_copy, edge[1], edge[2])
    end
    rem_edge!(graph_copy, start, finish)

    visited = Set{Int}()
    current_path = Vector{Int}()
    simple_paths = Vector{Vector{Int}}()
    get_simple_paths_dfs!(visited, start, finish, max_length, graph_copy,
                                current_path, simple_paths)
    path_bools = Vector{BoolPath}()
    for path in simple_paths
        bools::Vector{Bool} = [false for _ in 1:length(path)-1]
        for i in 1 : length(path)-1
            src_node = path[i]
            dst_node = path[i+1]
            if dst_node in outneighbors(query_graph, src_node)
                bools[i] = true # out edge
            else
                bools[i] = false # in edge
            end
        end
        push!(path_bools, bools)
    end

    includes_1_length_path = start in outneighbors(query_graph, finish)
    if includes_1_length_path
        push!(path_bools, Bool[false])
    end
    return path_bools
end

"""
Finds only the shortest path from the start to destination node in the Query Graph using
the A* search algorithm.
Returns a DiGraph representing the shortest path.
# Arguments
- start::Int - the starting node of the path in the Query Graph.
- finish::Int - the destination node of the path in the Query Graph.
- query::QueryGraph - the Query Graph to search for the shortest path.
"""
function get_matching_graph(start::Int, finish::Int, query::QueryGraph)
    # convert the graph to be undirected
    graph_copy = Graph(copy(query.graph))
    rem_edge!(graph_copy, start, finish)
    # get a path from the start to finish node
    edges = a_star(graph_copy, start, finish)
    # currently assumes the edges go in order of the path,
    # will have to debug through and see if it works
    new_graph = DiGraph(length(edges) + 1)
    current_start = 0
    for edge in edges
        current_start += 1
        if src(edge) in outneighbors(query.graph, dst(edge))
            # this is a backwards edge
            add_edge!(new_graph, current_start + 1, current_start)
        else
            # this is a forwards edge
            add_edge!(new_graph, current_start, current_start + 1)
        end
    end
    return new_graph
end

"""
For all partial paths, scales down their corresponding weights according the likelihood of a cycle closing with an edge from the end to the beginning of the path,
if there are any remaining unprocessed edges. This function is intended for usage once we finish traversing the query's spanning tree.
# Arguments
- query::QueryGraph - the Query Graph we are processing.
- summary::ColorSummary{DS} - the lifted summary of the Data Graph used for estimation.
- partial_paths::Array{Color} - the list of partial paths describing our current traversal through the Query Graph.
- partial_weights::Vector{W} - the list of weights corresponding to each partial path.
- current_query_nodes::Vector{Int} - the list of Query Nodes, describing the order we traverse through the Query Graph.
- visited_query_edges::Vector{Tuple{Int,Int}} - the list of edges (start, finish) that have already been processed.
- usingStoredStats::Bool - whether or not to use stored cycle statistics or assume independence during the estimation.
- only_shortest_path_cycle::Bool - whether to treat the closing edge as closing the shortest path (true) or as closing multiple paths of varying lengths (false)
"""
function handle_extra_edges!(query::QueryGraph, summary::ColorSummary{DS}, partial_paths::Array{Color}, partial_weights::Vector{W},
                                current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int,Int}}, usingStoredStats::Bool,
                                only_shortest_path_cycle::Bool) where DS where W
    # To account for cyclic queries, we check whether there are any remaining edges that have not
    # been processed. If so, we set the lower bound to 0, reduce the average estimate accordingly, and leave
    # the upper bound unchanged.
    remaining_edges::Vector{Tuple{Int, Int}} = []
    for edge in edges(query.graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges) &&
                 (src(edge) in current_query_nodes && dst(edge) in current_query_nodes)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end
    # scale down the average if there are remaining non-tree-edges
    for edge in remaining_edges
        push!(visited_query_edges, edge)
        parent_node_idx::Int = only(indexin(edge[1], current_query_nodes))
        child_node_idx::Int = only(indexin(edge[2], current_query_nodes))
        edge_label::Int = only(query.edge_labels[(edge[1],edge[2])])
        child_label::Int = only(query.vertex_labels[edge[2]])
        all_path_bools::Vector{BoolPath} = []
        if only_shortest_path_cycle
            all_path_bools = [convert_path_graph_to_bools(get_matching_graph(edge[1], edge[2], query))]
        else
            all_path_bools = get_all_simple_path_bools(edge[1], edge[2], summary.max_cycle_size, query.graph, visited_query_edges)
        end

        default_colors::StartEndColorPair = (-1, -1)
        path_counts = counter(all_path_bools)

        default_no_edge_probability::Float64 = 1.0
        default_no_edge_probabilities::Dict{BoolPath, Float64} = Dict()
        for (path_bools, path_count) in path_counts
            probability_no_edge = 1.0
            default_cycle_description = CyclePathAndColors(path_bools, default_colors)
            path_length = length(path_bools)
            if haskey(summary.cycle_probabilities, default_cycle_description)
                probability_no_edge *= (1.0 - summary.cycle_probabilities[default_cycle_description])^path_count
            elseif haskey(summary.cycle_length_probabilities, path_length)
                probability_no_edge *= (1.0 - summary.cycle_length_probabilities[path_length])^path_count
            else
                probability_no_edge *= (1.0 - get_independent_cycle_likelihood(summary))^path_count
            end

            default_no_edge_probability *= probability_no_edge
            default_no_edge_probabilities[path_bools] = probability_no_edge
        end


        edge_deg::Dict{Int, Dict{Int, DS}} = Dict()
        if haskey(summary.edge_deg, edge_label) && haskey(summary.edge_deg[edge_label], child_label)
            edge_deg = summary.edge_deg[edge_label][child_label]
        end
        for i  in 1:size(partial_paths)[2]
            parent_color::Color = partial_paths[parent_node_idx, i]
            child_color::Color =  partial_paths[child_node_idx, i]
            current_colors::StartEndColorPair = (parent_color, child_color)

            # We don't have to check data label because these nodes are already in the
            # partial path, so we have already ensured that the colors are appropriate
            probability_no_edge::Float64 = 1.0
            if (haskey(edge_deg, parent_color) && haskey(edge_deg[parent_color], child_color))
                if usingStoredStats && length(path_counts) > 0
                    for (path_bools, path_count) in path_counts
                        current_cycle_description = CyclePathAndColors(path_bools, current_colors)
                        if haskey(summary.cycle_probabilities, current_cycle_description)
                            probability_no_edge *= (1.0 - summary.cycle_probabilities[current_cycle_description])^path_count
                        else
                            probability_no_edge *= default_no_edge_probabilities[path_bools]
                        end
                    end
                else
                    probability_no_edge *= 1.0 - get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary)
                end
            end
            probability_no_edge *= 1.0 - summary.total_added_edges/summary.total_nodes^2
            # probability_no_edge *= ((summary.total_nodes^2-1)/(summary.total_nodes^2))^summary.total_added_edges
            partial_weights[i] = scale_coloring(partial_weights[i], (1.0 - probability_no_edge))
        end
    end
end

"""
Sums over all nodes that have already been processed when searching the Query Graph. Returns a condensed list of partial paths and their
corresponding weights.
# Arguments
- query::QueryGraph - the Query Graph that is currently being processed.
- partial_paths::Matrix{Color} - the matrix of partial paths describing the current traversal over the lifted color graph.
- partial_weights::Vector{W} - the list of weights for each partial path.
- current_query_nodes - the list of query nodes, describing the order we process the query graph.
- visited_query_edges::Vector{Tuple{Int, Int}} - the list of query edges that have already been processed.
"""
function sum_over_finished_query_nodes(query::QueryGraph, partial_paths::Matrix{Color}, partial_weights::Vector{W},
                                            current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int, Int}}) where W
    new_partial_paths, new_partial_weights = partial_paths, partial_weights
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
            new_partial_paths, new_partial_weights = sum_over_node(new_partial_paths, new_partial_weights, current_query_nodes, node)
        end
    end
    return new_partial_paths, new_partial_weights
end

"""
Estimates the cardinality of the Query Graph in the overall Data Graph using a lifted summary.
Finds all partial paths representing different color assignments for each query node then sums the weights
of all the partial paths to calculate the overall count. Returns a Float estimate.
# Arguments
- query::QueryGraph - the Query Graph to estimate the cardinality in the overall Data Graph.
- summary::ColorSummary{DS} - the summary storing statistics about the Data Graph after applying a graph coloring.
- max_partial_paths::Union{Nothing, Int} - the maximum number of partial paths to allow during estimation (enabling sampling).
- use_partial_sums::Bool - whether or not to sum over nodes during estimation.
- verbose::Bool - whether or not to output status messages during estimation.
- usingStoredStats::Bool - whether or not to use summary cycle statistics instead of assuming independence for cycle-closing likelihoods.
- include_cycles::Bool - whether or not to include cycles in the estimation process.
- sampling_strategy::SAMPLING_STRATEGY - the strategy to use for sampling partial paths.
- only_shortest_path_cycle::Bool - whether or not to treat cycle-closing edges as single-path closing events instead of closing multiple paths.
- timeout::Float64 - the maximum time to spend on estimation before returning a timeout error.
"""
function get_cardinality_bounds(query::QueryGraph, summary::ColorSummary{DS}; max_partial_paths::Union{Nothing, Int} = nothing,
                                use_partial_sums::Bool = true, verbose::Bool = false, usingStoredStats::Bool = false,
                                include_cycles::Bool = true, sampling_strategy::SAMPLING_STRATEGY=weighted,
                                only_shortest_path_cycle::Bool=false, timeout::Float64 = Inf) where DS
    start_time = time()
    W = stat_type_to_accumulator(DS)
    node_order::Vector{Int} = get_min_width_node_order(query.graph) #spanning tree to cut out cycles
    if verbose
        println("Node Order:", node_order)
    end
    # Because the label is implied by the color -> query_graph_vertex mapping stored in current_query_nodes,
    # we don't have to keep the label in the partial paths object.
    num_colors = summary.num_colors
    partial_paths = zeros(Color, 1, num_colors) # each tuple contains a pairing of color paths -> bounds
    partial_weights = fill(W(0), num_colors)
    visited_query_edges::Vector{Tuple{Int,Int}} = []
    current_query_nodes::Vector{Int} = []

    old_node = popfirst!(node_order)
    parent_label = only(query.vertex_labels[old_node])
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
                if data_label in summary.color_filters[color] || color in summary.color_full
                    data_label_is_in_color = true
                end
            end
            if (data_label_is_in_color)
                partial_paths[1, color] = color
                partial_weights[color] = W(summary.color_label_cardinality[color][parent_label])
            end
        end
    end
    new_node = old_node
    while length(node_order) > 0
        num_current_paths = size(partial_paths)[2]
        num_current_paths * num_colors > 10^6 && return get_default_count(DS) # If the requested memory is too large, return a timeout
        time() - start_time > timeout && return get_default_count(DS)
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        if use_partial_sums
            partial_paths, partial_weights = sum_over_finished_query_nodes(query, partial_paths, partial_weights, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end
        # Get the next child from the node order.
        new_node = popfirst!(node_order)
        parent_idx = 0
        out_edge = false
        for neighbor in all_neighbors(query.graph, new_node)
            if neighbor in current_query_nodes
                old_node = neighbor
                parent_idx::Int = only(indexin(neighbor, current_query_nodes))
                if old_node in inneighbors(query.graph, new_node)
                    out_edge = true
                end
                break
            end
        end
        # Push the current edge and nodes to the visited lists.
        if out_edge
            push!(visited_query_edges, (old_node, new_node))
        else
            push!(visited_query_edges, (new_node, old_node))
        end
        push!(current_query_nodes, new_node)

        # Get the appropriate labels, the query only uses one label per vertex/node.
        edge_label = 0
        if out_edge
            edge_label = only(query.edge_labels[(old_node,new_node)])
        else
            edge_label =  only(query.edge_labels[(new_node,old_node)])
        end

        new_label = only(query.vertex_labels[new_node])
        new_data_labels = get_data_label(query, new_node)
        num_current_paths = size(partial_paths)[2]
        num_current_paths * num_colors > 10^6 && return get_default_count(DS) # If the requested memory is too large, return a timeout.
        new_partial_paths = zeros(Color, length(current_query_nodes),  num_current_paths * num_colors)
        new_partial_weights = fill(W(0), num_current_paths * num_colors)
        # Update the partial paths using the parent-child combo that comes next from the query.
        edge_deg::Dict{Color, Dict{Color, DS}} = Dict()
        if haskey(summary.edge_deg, edge_label) &&
                        haskey(summary.edge_deg[edge_label], new_label)
            edge_deg = summary.edge_deg[edge_label][new_label]
        end
        new_path_idx = 1
        for i in 1:num_current_paths
            time() - start_time > timeout && return get_default_count(DS)
            old_color = partial_paths[parent_idx, i]
            # Account for colors with no outgoing children.
            if haskey(edge_deg, old_color)
                for new_color in keys(edge_deg[old_color])
                    time() - start_time > timeout && return get_default_count(DS)
                    # revamp the logic to use a set of labels rather than just one
                    # check if the data label(s) are in the color
                    data_label_in_color = false
                    for data_label in new_data_labels
                        if data_label == -1
                            data_label_in_color = true
                            continue
                        end
                        if data_label in summary.color_filters[new_color] || new_color in summary.color_full
                            data_label_in_color = true
                        end
                    end
                    if !data_label_in_color
                        continue
                    end
                    for j in 1:length(current_query_nodes)-1
                        new_partial_paths[j, new_path_idx] = partial_paths[j, i]
                    end
                    new_partial_paths[length(current_query_nodes), new_path_idx] = new_color
                    degree_stats::DS = edge_deg[old_color][new_color]
                    new_partial_weights[new_path_idx] = extend_coloring(partial_weights[i], degree_stats, out_edge)
                    if !(length(new_data_labels) == 1 && new_data_labels[1] == -1)
                        new_partial_weights[new_path_idx] = scale_coloring(new_partial_weights[new_path_idx], 1.0/summary.color_label_cardinality[new_color][new_label])
                    end
                    isinf(get_count(new_partial_weights[new_path_idx])) && return Inf
                    new_path_idx += 1
                end
            end
        end
        time() - start_time > timeout && return get_default_count(DS)

        partial_paths = new_partial_paths[:, 1:new_path_idx-1]
        partial_weights = new_partial_weights[1:new_path_idx-1]
        if !(isnothing(max_partial_paths)) && (size(partial_paths)[2] > max_partial_paths)
            partial_paths, partial_weights = sample_paths(partial_paths, partial_weights, max_partial_paths, sampling_strategy)
        end

        if (include_cycles)
            handle_extra_edges!(query, summary, partial_paths, partial_weights, current_query_nodes, visited_query_edges, usingStoredStats, only_shortest_path_cycle)
        end
    end

    # Sum over the calculated partial paths to get the final bounds.
    final_estimate = 0.0
    for w in partial_weights
        final_estimate += get_count(w)
    end
    return final_estimate
end
