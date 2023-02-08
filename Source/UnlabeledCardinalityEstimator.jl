
# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
# It currently only handles query graphs without labels.
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
using QuasiStableColors

struct ColorSummary
    color_cardinality::Dict{Int32, Int32}
    edge_cardinality::Dict{Int32, Dict{Int32, Float64}}
    edge_min_deg::Dict{Int32, Dict{Int32, Float64}}
    edge_avg_deg::Dict{Int32, Dict{Int32, Float64}}
    edge_max_deg::Dict{Int32, Dict{Int32, Float64}}
end

@noinline function generate_color_summary(g::DiGraph, numColors::Int)
    color_cardinality = counter(Int)
    undirected_graph = Graph(g)
    C = q_color(undirected_graph, n_colors=numColors)
    color_hash::Dict{Int, Int32} = Dict()
    for (color, nodes) in enumerate(C)
        for x in nodes
            color_hash[x] = color
            inc!(color_cardinality, color)
        end
    end

    color_to_color_counter::Dict{Int32, Dict{Int32, Any}} = Dict()
    for x in vertices(g)
        c1 = color_hash[x]
        for y in outneighbors(g,x) 
            c2 = color_hash[y]
            if !haskey(color_to_color_counter, c1)
                color_to_color_counter[c1] = Dict()
            end
            if !haskey(color_to_color_counter[c1], c2)
                color_to_color_counter[c1][c2] = counter(Int32)
            end
            inc!(color_to_color_counter[c1][c2], x)
        end
    end

    edge_cardinality::Dict{Int32, Dict{Int32, Float64}} = Dict()
    edge_min_deg::Dict{Int32, Dict{Int32, Float64}} = Dict()
    edge_avg_deg::Dict{Int32, Dict{Int32, Float64}} = Dict()
    edge_max_deg::Dict{Int32, Dict{Int32, Float64}} = Dict()
    for c1 in keys(color_to_color_counter)
        edge_min_deg[c1] = Dict()
        edge_avg_deg[c1] = Dict()
        edge_max_deg[c1] = Dict()
        edge_cardinality[c1] = Dict()
        for c2 in keys(color_to_color_counter[c1])
            edge_min_deg[c1][c2] =  nv(g) 
            edge_max_deg[c1][c2] =  0 
            edge_cardinality[c1][c2] = 0
            for v in values(color_to_color_counter[c1][c2])
                edge_min_deg[c1][c2] =  min(v, edge_min_deg[c1][c2]) 
                edge_max_deg[c1][c2] = max(v, edge_max_deg[c1][c2])
                edge_cardinality[c1][c2] += v
            end

            # `color_cardinality` is used here to account for nodes in c1 without any edges to c2.
            edge_avg_deg[c1][c2] = edge_cardinality[c1][c2]/color_cardinality[c1]
            if length(values(color_to_color_counter[c1][c2])) < color_cardinality[c1]
                edge_min_deg[c1][c2] =  0 
            end
        end
    end

    return ColorSummary(color_cardinality, edge_cardinality, edge_min_deg, edge_avg_deg, edge_max_deg)
end 

# The following two functions sum over all paths which have the same color assigned to a particular node in the query graph.
# Equivalently, they perform a groupby on all other nodes of the query graph. The goal of this is to prevent
# an exponential growth in the number of paths through the lifted color graph. However, we can only remove query nodes whose
# edges have already been processed.
function sum_over_node!(partial_paths, current_query_nodes, node_to_remove)
    nodeIdx = 1
    for node in current_query_nodes
        if node == node_to_remove
            break
        end 
        nodeIdx += 1
    end
    new_partial_paths::Dict{Array{Int}, Union{Array{Float64}, Int}} = Dict()
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

function sum_over_finished_query_nodes!(query_graph, partial_paths, current_query_nodes, visited_query_edges)
    prev_query_nodes = copy(current_query_nodes)
    for node in prev_query_nodes
        has_living_edges = false
        for neighbor in outneighbors(query_graph, node)
            if !((node, neighbor) in visited_query_edges)
                has_living_edges = true
            end
        end
        for neighbor in inneighbors(query_graph, node)
            if !((neighbor, node) in visited_query_edges)
                has_living_edges = true
            end
        end
        if ! has_living_edges
            sum_over_node!(partial_paths, current_query_nodes, node)
        end
    end
end

function get_cardinality_bounds_given_starting_node(query_graph::DiGraph, summary::ColorSummary, 
                                                    starting_node::Int; use_partial_sums = true, verbose = false)
    node_order = topological_sort_by_dfs(bfs_tree(query_graph, starting_node))
    partial_paths::Array{Tuple{Array{Int}, Array{Float64}}} = []
    visited_query_edges = []
    current_query_nodes = []
    parent_node = popfirst!(node_order)
    child_node = popfirst!(node_order)
    push!(visited_query_edges, (parent_node, child_node))
    push!(current_query_nodes, parent_node)
    push!(current_query_nodes, child_node)
    for c1 in keys(summary.edge_cardinality)
        for c2 in keys(summary.edge_cardinality[c1])
            push!(partial_paths, ([c1,c2], [summary.edge_cardinality[c1][c2],
                                     summary.edge_cardinality[c1][c2],
                                     summary.edge_cardinality[c1][c2]]))
        end
    end

    while length(node_order) > 0
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        if use_partial_sums
            sum_over_finished_query_nodes!(query_graph, partial_paths, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end

        child_node = popfirst!(node_order)
        parent_idx = 0
        for neighbor in inneighbors(query_graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
            end
        end
        push!(visited_query_edges, (parent_node, child_node))
        push!(current_query_nodes, child_node)

        new_partial_paths::Array{Tuple{Array{Int}, Array{Float64}}} = []
        for path_and_bounds in partial_paths
            path = path_and_bounds[1]
            running_bounds = path_and_bounds[2]
            parent_color = only(path[parent_idx])
            # account for colors with no outgoing children
            if (!haskey(summary.edge_avg_deg, parent_color))
                continue
            end
            for child_color in keys(summary.edge_avg_deg[parent_color])
                new_path = copy(path)
                push!(new_path, child_color)
                new_bounds = [running_bounds[1]*summary.edge_min_deg[parent_color][child_color],
                              running_bounds[2]*summary.edge_avg_deg[parent_color][child_color],
                              running_bounds[3]*summary.edge_max_deg[parent_color][child_color],
                ]
                push!(new_partial_paths, (new_path, new_bounds))
            end
        end
        partial_paths = new_partial_paths
    end

    # To account for cyclic queries, we check whether there are any remaining edges that have not
    # been processed. If so, we set the lower bound to 0, reduce the average estimate accordingly, and leave
    # the upper bound unchanged.
    remaining_edges = []
    for edge in edges(query_graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end

    final_bounds = [0,0,0]
    for path_and_bounds in partial_paths
        path = path_and_bounds[1]
        bounds = path_and_bounds[2]
        lower = only(bounds[1])
        if length(remaining_edges) > 0
            lower = 0
        end
        average = only(bounds[2])
        for edge in remaining_edges
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_color = only(path[parent_node_idx])
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_color = only(path[child_node_idx])
            probability_of_edge = 0
            if haskey(summary.edge_avg_deg, parent_color) && haskey(summary.edge_avg_deg[parent_color], child_color)
                probability_of_edge = summary.edge_avg_deg[parent_color][child_color]/summary.color_cardinality[child_color]
            end
            average *= probability_of_edge
        end
        
        upper = only(bounds[3])
        final_bounds = final_bounds .+ [lower, average, upper]
    end 
    return final_bounds
end

function get_cardinality_bounds(query_graph::DiGraph, summary::ColorSummary; 
                                use_partial_sums = true, try_all_starting_nodes=true, verbose = false)
    root_nodes = []
    for node in vertices(query_graph)
        if is_connected(bfs_tree(query_graph, node))
            push!(root_nodes, node)
        end
    end
    if try_all_starting_nodes
        final_bounds = [0, 0, Inf]
        for node in root_nodes
            cur_bounds = get_cardinality_bounds_given_starting_node(query_graph, summary, node;
                                                                    use_partial_sums=use_partial_sums, verbose=verbose)
            final_bounds[1] = max(final_bounds[1], cur_bounds[1])
            final_bounds[2] += cur_bounds[2]
            final_bounds[3] = min(final_bounds[3], cur_bounds[3])
        end
        final_bounds[2] /= length(root_nodes)
        return final_bounds
    end
    return get_cardinality_bounds_given_starting_node(query_graph, summary, root_nodes[1];
                                                      use_partial_sums=use_partial_sums, verbose=verbose)
end

# We use the same general structure to calculate the exact size of the query by finding all paths
# on the original data graph and giving each path a weight of 1.
function get_exact_size(query_graph::DiGraph, data_graph::DiGraph; use_partial_sums = true, verbose=false)
    node_order = topological_sort_by_dfs(bfs_tree(query_graph, vertices(query_graph)[1]))
    partial_paths::Array{Tuple{Array{Int}, Int}} = []
    visited_query_edges = []
    current_query_nodes = []
    parent_node = popfirst!(node_order)
    child_node = popfirst!(node_order)
    push!(visited_query_edges, (parent_node, child_node))
    push!(current_query_nodes, parent_node)
    push!(current_query_nodes, child_node)
    for c1 in vertices(data_graph)
        for c2 in outneighbors(data_graph, c1)
            push!(partial_paths, ([c1,c2], 1))
        end
    end

    while length(node_order) > 0
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        if use_partial_sums
            sum_over_finished_query_nodes!(query_graph, partial_paths, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end
        
        child_node = popfirst!(node_order)
        parent_idx = 0
        for neighbor in inneighbors(query_graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
            end
        end
        
        push!(current_query_nodes, child_node)
        push!(visited_query_edges, (parent_node, child_node))
        new_partial_paths::Array{Tuple{Array{Int}, Int}} = []
        for path_and_weight in partial_paths
            path = path_and_weight[1]
            weight = path_and_weight[2]
            parent_node = only(path[parent_idx])
            for child_node in outneighbors(data_graph, parent_node)
                new_path = copy(path)
                push!(new_path, child_node)
                push!(new_partial_paths, (new_path, weight))
            end
        end
        partial_paths = new_partial_paths
    end
    remaining_edges = []
    for edge in edges(query_graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end
    
    final_bounds = 0
    for path_and_weight in partial_paths 
        path = path_and_weight[1]
        weight = path_and_weight[2]
        satisfies_cycles = true
        for edge in remaining_edges
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_data_node = only(path[parent_node_idx])
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_data_node = only(path[child_node_idx])
            if !(child_data_node in inneighbors(data_graph, parent_data_node))
                satisfies_cycles = false
            end
        end 
        if satisfies_cycles
            final_bounds = final_bounds + weight
        end
    end 
    return final_bounds
end
