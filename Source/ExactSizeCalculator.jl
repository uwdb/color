# This file contains a prototype implementation of exact sub-graph counting.
include("PropertyGraph.jl")
using Graphs
using QuasiStableColors
using Probably
using StatsBase


function sum_over_node_exact!(partial_paths::Vector{Tuple{Vector{Int},  Int}}, current_query_nodes, node_to_remove)
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


function sum_over_finished_query_nodes_exact!(query::QueryGraph, partial_paths::Vector{Tuple{Vector{Int},  Int}}, current_query_nodes, visited_query_edges)
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
            sum_over_node_exact!(partial_paths, current_query_nodes, node)
        end
    end
end

function handle_extra_edges_exact!(query::QueryGraph, data::DataGraph, partial_paths::Vector{Tuple{Vector{Int}, Int}}, current_query_nodes, visited_query_edges)
    new_partial_paths::Vector{Tuple{Vector{Int}, Int}} = []
    remaining_edges = []
    for edge in edges(query.graph)
        # since the edge's nodes are already processed, we don't have to check 
        if ! ((src(edge), dst(edge)) in visited_query_edges) &&
                 (src(edge) in current_query_nodes && dst(edge) in current_query_nodes)
            push!(remaining_edges, (src(edge), dst(edge)))
            push!(visited_query_edges, (src(edge), dst(edge)))
        end
    end
    for path_and_weight in partial_paths 
        path = path_and_weight[1]
        weight = path_and_weight[2]
        satisfies_cycles = true
        for edge in remaining_edges
            # Only count the cycle as satisfied if this remaining edge's label matches the query graph's edge label.

            # Get the parent node from the list of current query nodes.
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_data_node = only(path[parent_node_idx])
            # Get the child node from the list of current query nodes.
            new_node_idx = indexin(edge[2], current_query_nodes)
            child_data_node = only(path[new_node_idx])

            # Check if the edge label exists, if it doesn't then we can break here.
            # Don't need to check parent node because we got the parent node from the data graph,
            # but we do need to check if there is an edge connection to the child.
            if (!haskey(data.edge_labels, (parent_data_node, child_data_node)))
                satisfies_cycles = false;
                break;
            end
            data_edge_labels = data.edge_labels[(parent_data_node,child_data_node)]
            data_child_vertex_labels = data.vertex_labels[child_data_node]
            query_edge_label = only(query.edge_labels[(edge[1],edge[2])])
            query_child_vertex_label = only(query.vertex_labels[edge[2]])
            if !((query_edge_label == -1 || in(query_edge_label, data_edge_labels)) && 
                    (query_child_vertex_label == -1 || in(query_child_vertex_label, data_child_vertex_labels)))
                satisfies_cycles = false
                break
            end
        end 
        if satisfies_cycles
            push!(new_partial_paths, (path, weight))
        end
    end 

    empty!(partial_paths)
    copy!(partial_paths, new_partial_paths)
    empty!(new_partial_paths)
end

# We use the same general structure to calculate the exact size of the query by finding all paths
# on the original data graph and giving each path a weight of 1. 
function get_exact_size(query::QueryGraph, data::DataGraph; use_partial_sums = true, verbose=false, starting_nodes = nothing, ending_nodes = nothing)
    node_order = get_min_width_node_order(query.graph)
    partial_paths::Vector{Tuple{Vector{Int}, Int}} = []
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
    push!(current_query_nodes, old_node)
    if starting_nodes === nothing
        starting_nodes = vertices(data.graph)
    end
    for node in starting_nodes
        # if the id labels don't match, then don't initialize with this node
        query_data_labels = get_data_label(query, new_node)
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
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges)
        if use_partial_sums
            sum_over_finished_query_nodes_exact!(query, partial_paths, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end

        new_node = popfirst!(node_order)
        parent_idx = 0
        outEdge = false
        for neighbor in all_neighbors(query.graph, new_node)
            if neighbor in current_query_nodes
                old_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
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
        query_child_label = query.vertex_labels[new_node][1]
        query_child_id_labels = query.vertex_id_labels[new_node]
        
        push!(current_query_nodes, new_node)
        new_partial_paths::Vector{Tuple{Vector{Int}, Int}} = []
        for path_and_weight in partial_paths
            path = path_and_weight[1]
            weight = path_and_weight[2]
            old_node = only(path[parent_idx])
            if outEdge
                for data_new_node in outneighbors(data.graph, old_node)
                    new_weight = weight
                    # Only add a new partial path if the edge label and node label match our query.
                    data_edge_labels = data.edge_labels[(old_node,data_new_node)]
                    data_child_labels = data.vertex_labels[data_new_node]
                    data_child_id_label = get_data_label(data, data_new_node)
                    if (query_child_id_labels != [-1] && length(intersect(query_child_id_labels,data_child_id_label)) == 0)
                        continue
                    end
                    if (query_child_label == -1  || in(query_child_label, data_child_labels)) && 
                        (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                        new_path = copy(path)
                        push!(new_path, data_new_node)
                        push!(new_partial_paths, (new_path, new_weight))
                    end
                end
            else
                for data_new_node in inneighbors(data.graph, old_node)
                    if (!(ending_nodes === nothing) && !(data_new_node in ending_nodes))
                    end
                    new_weight = weight
                    # Only add a new partial path if the edge label and node label match our query.
                    data_edge_labels = data.edge_labels[(data_new_node,old_node)]
                    data_child_labels = data.vertex_labels[data_new_node]
                    data_child_id_label = get_data_label(data, data_new_node)
                    if (query_child_id_labels != [-1])
                        if (length(intersect(query_child_id_labels, data_child_id_label)) == 0)
                            continue
                        end
                    end
                    if (query_child_label == -1  || in(query_child_label, data_child_labels)) && 
                        (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                        new_path = copy(path)
                        push!(new_path, data_new_node)
                        push!(new_partial_paths, (new_path, new_weight))
                    end
                end
            end
        end
        partial_paths = new_partial_paths
    end
    handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges)
    
    exact_size = 0
    for path_and_weight in partial_paths 
        exact_size +=  path_and_weight[2]
    end 
    return exact_size
end
