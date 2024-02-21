# This file contains a prototype implementation of exact sub-graph counting.

"""
Sums out a specific node from a list of paths, returning the condensed list of paths without the node. 
Used as an optimization for speed and storage during estimation time.
# Arguments
- partial_paths::Vector{Tuple{Vector{NodeId},  Int}} - a list of partial paths to condense, consisting of tuples of nodes in the path and the overall weight.
- current_query_nodes - the list of query nodes we are currently processing.
- node_to_remove - the node to sum over (and remove when condensing the paths).
- timeout - the maximum time to spend on this operation before returning an error.
- start_time - the time that the operation began.
"""
function sum_over_node_exact!(partial_paths::Vector{Tuple{Vector{NodeId},  Int}},
                                 current_query_nodes, node_to_remove,
                                 timeout, start_time)
    nodeIdx = 1
    for node in current_query_nodes
        if node == node_to_remove
            break
        end
        nodeIdx += 1
    end
    new_partial_paths::Dict{Vector{NodeId}, Union{Vector{Float64}, Int}} = Dict()
    for path_and_bounds in partial_paths
        if timeout > 0 && time() - start_time > timeout
            println("Timeout Reached")
            return
        end
        path = path_and_bounds[1]
        bounds = path_and_bounds[2]
        new_path = [path[begin:nodeIdx-1]..., path[nodeIdx+1:end]...]
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

"""
Sums over all nodes with no remaining connections in the query graph, condensing the list of paths.
# Arguments
- query::QueryGraph - the query graph that the partial paths are traversing.
- partial_paths::Vector{Tuple{Vector{NodeId}, Int}} - the list of partial paths describing the current estimation through the query graph.
- current_query_nodes - the list of query nodes, describing the order we are traversing through the query graph.
- visited_query_edges - the list of edges that have already been processed.
- timeout - the maximum time to spend on this operation before returning a timeout error.
- start_time - the time that the operation began.
- nodes_to_not_sum - specific nodes to avoid summing out, if any. 
"""
function sum_over_finished_query_nodes_exact!(query::QueryGraph, partial_paths::Vector{Tuple{Vector{NodeId},  Int}},
                                                 current_query_nodes, visited_query_edges,
                                                 timeout, start_time; nodes_to_not_sum = [])
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
        if !(has_living_edges) && !(node in nodes_to_not_sum)
            sum_over_node_exact!(partial_paths, current_query_nodes, node, timeout, start_time)
        end
    end
end

"""
Scales down the partial path results based on the remaining edges not included in the tree traversal (i.e. those that close a cycle).
# Arguments
- query::QueryGraph - the query graph to process.
- data::DataGraph - the data graph to find query graph matches in.
- partial_paths::Vector{Tuple{Vector{NodeId}, Int}} - the current list of paths traversing through the data graph.
- current_query_nodes - a list of nodes describing the order to process the nodes in the query graph.
- visited_query_edges - a list of the query edges which have already been processed.
- timeout - the maximum time to spend on this operation before returning a timeout error.
- start_time - the time that this function began.
"""
function handle_extra_edges_exact!(query::QueryGraph, data::DataGraph,
                                    partial_paths::Vector{Tuple{Vector{NodeId}, Int}},
                                    current_query_nodes, visited_query_edges,
                                    timeout, start_time)
    remaining_edges::Vector{Tuple{NodeId, NodeId}} = []
    for edge in edges(query.graph)
        # since the edge's nodes are already processed, we don't have to check
        if ! ((src(edge), dst(edge)) in visited_query_edges) &&
                 (src(edge) in current_query_nodes && dst(edge) in current_query_nodes)
            push!(remaining_edges, (src(edge), dst(edge)))
            push!(visited_query_edges, (src(edge), dst(edge)))
        end
    end
    if length(remaining_edges) == 0
        return
    end
    for i in eachindex(partial_paths)
        if timeout > 0 && time() - start_time > timeout
            println("Timeout Reached")
            return
        end
        path = partial_paths[i][1]
        weight = partial_paths[i][2]
        for edge in remaining_edges
            # Only count the cycle as satisfied if this remaining edge's label matches the query graph's edge label.

            # Get the parent node from the list of current query nodes.
            parent_node_idx = findfirst(x -> x == edge[1], current_query_nodes)
            parent_data_node = path[parent_node_idx]
            # Get the child node from the list of current query nodes.
            new_node_idx = findfirst(x -> x == edge[2], current_query_nodes)
            child_data_node = path[new_node_idx]

            # Check if the edge label exists, if it doesn't then we can break here.
            # Don't need to check parent node because we got the parent node from the data graph,
            # but we do need to check if there is an edge connection to the child.
            if (!haskey(data.edge_labels, (parent_data_node, child_data_node)))
                weight = 0
                break;
            end
            data_edge_labels = data.edge_labels[(parent_data_node,child_data_node)]
            data_child_vertex_labels = data.vertex_labels[child_data_node]
            query_edge_label = only(query.edge_labels[(edge[1],edge[2])])
            query_child_vertex_label = only(query.vertex_labels[edge[2]])
            if !((query_edge_label == -1 || in(query_edge_label, data_edge_labels)) &&
                    (query_child_vertex_label == -1 || in(query_child_vertex_label, data_child_vertex_labels)))
                weight = 0
                break
            end
        end
        partial_paths[i] = (path, weight)
    end
    filter!(x -> x[2] > 0, partial_paths)
end

"""
Returns a uniformly-sampled subset of the current partial paths then redistributes the weights of the unsampled paths
based on the magnitude of the remaining paths' weights.
# Arguments
- partial_paths - a list of partial paths describing the current traversal for the estimation.
- num_samples - the number of paths to keep.
"""
function sample_paths_exact(partial_paths, num_samples::Int)
    # partial_path[x] = (color path, bounds)

    # if we want to sample more paths than there are existing, then just return the original partial paths
    if (num_samples > length(partial_paths))
        return partial_paths
    end

    # sum up all of the bounds
    overall_bounds_sum = 0
    for path_and_bounds in partial_paths
        overall_bounds_sum = overall_bounds_sum + path_and_bounds[2]
    end

    # choose a sample of the paths
    path_samples::Vector{Tuple{Vector{NodeId}, Int}} = sample(partial_paths, num_samples; replace=false)

    # sum up the sampled bounds
    sampled_bounds_sum = 0
    for path_and_bounds in path_samples
        sampled_bounds_sum = sampled_bounds_sum + path_and_bounds[2]
    end
    # get the difference between the overall and sampled bound sum_over_finished_query_nodes
    bound_diff = overall_bounds_sum - sampled_bounds_sum

    # for each sampled path...
    for i in eachindex(path_samples)
        # figure out what fraction of the sampled bounds is in the current Bounds
        # higher bounds will have more weight redistributed to them
        bound_fractions = path_samples[i][2] / sampled_bounds_sum
        # use that fraction of the difference (i.e. the removed path weights) and add it to the partial path
        redistributed_weights = ceil(bound_fractions * bound_diff)
        path_samples[i] = (path_samples[i][1], path_samples[i][2] + redistributed_weights)
    end

    return path_samples
end

"""
Uses the same general structure as the estimation to calculate the exact size of the query by finding all paths
on the original data graph and giving each path a weight of 1. Returns a list of partial paths where each path is a match
of the query graph in the overall data graph.
# Arguments
- query::QueryGraph - the query graph to process.
- data::DataGraph - the data graph to find query graph matches in.
- use_partial_sums - whether or not to sum over finished nodes to optimize runtime and storage.
- verbose - whether or not to print status updates during the counting.
- max_partial_paths - the maximum number of partial paths to use during processing, enabling sampling during counting.
- timeout - the maximum time to spend on this operation before returning a timeout error.
- nodes_to_keep - specific nodes to avoid summing out, so the method can return specific paths.
"""
function get_subgraph_counts(query::QueryGraph, data::DataGraph; use_partial_sums = true, verbose=false,
                                                                max_partial_paths = nothing,
                                                                timeout = -1, nodes_to_keep = [])
    start_time = time()

    node_order = get_min_width_node_order(query.graph)
    partial_paths::Vector{Tuple{Vector{NodeId}, Int}} = []
    visited_query_edges = []
    current_query_nodes = []
    if verbose
        println("Node Order: ", node_order)
    end
    # Initialize partial_paths as 1-node paths with label matching the
    # initial query node's label.
    old_node = popfirst!(node_order)
    new_node = old_node
    parent_label = only(query.vertex_labels[old_node])
    query_data_labels = get_data_label(query, new_node)
    push!(current_query_nodes, old_node)
    for node in 1:nv(data.graph)
        # if the id labels don't match, then don't initialize with this node
        if (query_data_labels != [-1] && length(intersect(query_data_labels, get_data_label(data, node))) == 0)
            continue
        end
        node_labels = data.vertex_labels[node]
        # if the node labels don't match, then don't initialize with this node
        if (parent_label == -1) || (in(parent_label, node_labels))
            push!(partial_paths, ([node], 1))
        end
    end


    while length(node_order) > 0
        if use_partial_sums
            # pass nodes_to_keep here, just don't sum over it
            sum_over_finished_query_nodes_exact!(query, partial_paths, current_query_nodes, visited_query_edges, timeout, start_time, nodes_to_not_sum = nodes_to_keep)
        end

        if !(isnothing(max_partial_paths)) && (length(partial_paths) > max_partial_paths)
            partial_paths = sample_paths_exact(partial_paths, max_partial_paths)
        end

        new_node = popfirst!(node_order)
        parent_idx::Int = 0
        outEdge = false
        for neighbor in all_neighbors(query.graph, new_node)
            if neighbor in current_query_nodes
                old_node = neighbor
                parent_idx = findfirst(x -> x == neighbor, current_query_nodes)
                if neighbor in inneighbors(query.graph, new_node)
                    outEdge = true
                end
                break
            end
        end
        query_edge_label = 0
        if outEdge
            query_edge_label = only(query.edge_labels[(old_node,new_node)])
            push!(visited_query_edges, (old_node, new_node))
        else
            query_edge_label =  only(query.edge_labels[(new_node,old_node)])
            push!(visited_query_edges, (new_node, old_node))
        end
        query_child_label = only(query.vertex_labels[new_node])
        query_child_id_labels = query.vertex_id_labels[new_node]

        push!(current_query_nodes, new_node)
        new_partial_paths::Vector{Tuple{Vector{NodeId}, Int}} = []
        sizehint!(new_partial_paths, length(partial_paths))
        for path_and_weight in partial_paths
            if timeout > 0 && time() - start_time > timeout
                println("Timeout Reached")
                return -1
            end
            path::Vector{NodeId} = path_and_weight[1]
            weight = path_and_weight[2]
            old_node = path[parent_idx]
            wildcard_label = [-1]
            if outEdge
                if !isnothing(max_partial_paths)
                    valid_neighbors = 0
                    for data_new_node in outneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(old_node, data_new_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                                continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            valid_neighbors += 1
                        end
                    end
                    valid_neighbors == 0 && continue
                    sampled_node = rand(1:valid_neighbors)
                    node_count = 0
                    for data_new_node in outneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(old_node, data_new_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                            continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            node_count += 1
                            if node_count == sampled_node
                                new_path::Vector{NodeId} = copy(path)
                                push!(new_path, data_new_node)
                                push!(new_partial_paths, (new_path, weight*valid_neighbors))
                                break
                            end
                        end
                    end
                else
                    for data_new_node in outneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(old_node, data_new_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                            continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            new_path::Vector{NodeId} = copy(path)
                            push!(new_path, data_new_node)
                            push!(new_partial_paths, (new_path, weight))
                        end
                    end
                end
            else
                if !isnothing(max_partial_paths)
                    valid_neighbors = 0
                    for data_new_node in inneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(data_new_node,old_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                                continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            valid_neighbors += 1
                        end
                    end
                    valid_neighbors == 0 && continue
                    sampled_node = rand(1:valid_neighbors)
                    node_count = 0
                    for data_new_node in inneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(data_new_node,old_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                            continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            node_count += 1
                            if node_count == sampled_node
                                new_path::Vector{NodeId} = copy(path)
                                push!(new_path, data_new_node)
                                push!(new_partial_paths, (new_path, weight*valid_neighbors))
                                break
                            end
                        end
                    end
                else
                    for data_new_node in inneighbors(data.graph, old_node)
                        # Only add a new partial path if the edge label and node label match our query.
                        data_edge_labels = data.edge_labels[(data_new_node,old_node)]
                        data_child_labels = data.vertex_labels[data_new_node]
                        data_child_id_label = get_data_label(data, data_new_node)
                        if (query_child_id_labels != wildcard_label) && (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                            continue
                        end
                        if (query_child_label == -1  || in(query_child_label, data_child_labels)) &&
                            (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                            new_path::Vector{NodeId} = copy(path)
                            push!(new_path, data_new_node)
                            push!(new_partial_paths, (new_path, weight))
                        end
                    end
                end
            end
        end
        partial_paths = new_partial_paths

        if !(isnothing(max_partial_paths)) && (length(partial_paths) > max_partial_paths)
            partial_paths = sample_paths_exact(partial_paths, max_partial_paths)
        end

        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges, timeout, start_time)
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end
    end
    handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges, timeout, start_time)

    # have one final call to sum over sum_over_finished_query_nodes_exact
    sum_over_finished_query_nodes_exact!(query, partial_paths, current_query_nodes, visited_query_edges, timeout, time(), nodes_to_not_sum = nodes_to_keep)

    return partial_paths # if we sum over everything, this will be empty
end

"""
Calculates the exact size of the number of query graph appearances in the data graph by finding all partial paths of query matches and summing their counts.
# Arguments
- query::QueryGraph - the query graph to process.
- data::DataGraph - the data graph to find query graph matches in.
- use_partial_sums - whether or not to sum over finished nodes to optimize runtime and storage.
- verbose - whether or not to print status updates during the counting.
- max_partial_paths - the maximum number of partial paths to use during processing, enabling sampling during counting.
- timeout - the maximum time to spend on this operation before returning a timeout error.
- nodes_to_keep - specific nodes to avoid summing out, so the method can return specific paths.
"""
function get_exact_size(query::QueryGraph, data::DataGraph; use_partial_sums = true, verbose=false, max_partial_paths = nothing, timeout = -1, nodes_to_keep = [])
    partial_paths = get_subgraph_counts(query, data; use_partial_sums=use_partial_sums, verbose=verbose, max_partial_paths = max_partial_paths, timeout=timeout, nodes_to_keep = nodes_to_keep)
    if partial_paths  == -1
        return -1
    end
    exact_size = 0
    for path_and_weight in partial_paths
        exact_size += path_and_weight[2]
    end
    return exact_size
end
