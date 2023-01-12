# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
# It currently only handles acyclic query graphs without labels.
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
using QuasiStableColors

struct ColorSummary
    edge_cardinality::Dict{Int, Dict{Int, Float64}}
    edge_min_deg::Dict{Int, Dict{Int, Float64}}
    edge_avg_deg::Dict{Int, Dict{Int, Float64}}
    edge_max_deg::Dict{Int, Dict{Int, Float64}}
end

function generate_color_summary(g::Graph, numColors::Int)
    color_cardinality = counter(Int)
    C = q_color(g, n_colors=numColors)
    color_hash::Dict{Int, Int} = Dict()
    for (color, nodes) in enumerate(C)
        for x in nodes
            color_hash[x] = color
            inc!(color_cardinality, color)
        end
    end

    color_to_color_counter::Dict{Int, Dict{Int, counter}} = Dict()
    for x in vertices(g)
        c1 = color_hash[x]
        for y in all_neighbors(g,x) 
            c2 = color_hash[y]
            if !haskey(color_to_color_counter, c1)
                color_to_color_counter[c1] = Dict()
            end
            if !haskey(color_to_color_counter[c1], c2)
                color_to_color_counter[c1][c2] = counter(Int)
            end
            inc!(color_to_color_counter[c1][c2], x)
        end
    end

    edge_cardinality::Dict{Int, Dict{Int, Float64}} = Dict()
    edge_min_deg::Dict{Int, Dict{Int, Float64}} = Dict()
    edge_avg_deg::Dict{Int, Dict{Int, Float64}} = Dict()
    edge_max_deg::Dict{Int, Dict{Int, Float64}} = Dict()
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

    return ColorSummary(edge_cardinality, edge_min_deg, edge_avg_deg, edge_max_deg)
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
    new_partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
    for path in keys(partial_paths)
        new_path = copy(path)
        deleteat!(new_path, nodeIdx)
        if ! haskey(new_partial_paths, new_path)
            new_partial_paths[new_path] = partial_paths[path]
        else
            new_partial_paths[new_path] = new_partial_paths[new_path] .+ partial_paths[path]
        end
    end
    deleteat!(current_query_nodes, nodeIdx)
    empty!(partial_paths)
    for path in keys(new_partial_paths)
        partial_paths[path] = new_partial_paths[path]
    end
end

function sum_over_finished_query_nodes!(query_graph, partial_paths, current_query_nodes, visited_query_nodes)
    prev_query_nodes = copy(current_query_nodes)
    for node in prev_query_nodes
        has_living_neighbors = false
        for neighbor in all_neighbors(query_graph, node)
            if !(neighbor in visited_query_nodes)
                has_living_neighbors = true
            end
        end
        if !has_living_neighbors & use_partial_sums
            sum_over_node!(partial_paths, current_query_nodes, node)
        end
    end
end

function get_cardinality_bounds(query_graph::DiGraph, summary::ColorSummary; use_partial_sums = true, verbose = false)
    if is_cyclic(query_graph)
        println("Query Graph Must Be Acyclic")
        return -1
    end

    node_order = topological_sort_by_dfs(query_graph)
    partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
    visited_query_nodes = []
    current_query_nodes = []
    parent_node = popfirst!(node_order)
    child_node = popfirst!(node_order)
    push!(visited_query_nodes, parent_node)
    push!(visited_query_nodes, child_node)
    push!(current_query_nodes, parent_node)
    push!(current_query_nodes, child_node)
    for c1 in keys(summary.edge_cardinality)
        for c2 in keys(summary.edge_cardinality[c1])
            partial_paths[[c1,c2]] = [summary.edge_cardinality[c1][c2],
                                     summary.edge_cardinality[c1][c2],
                                     summary.edge_cardinality[c1][c2]]
        end
    end

    while length(node_order) > 0
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Nodes: ", visited_query_nodes)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        sum_over_finished_query_nodes!(query_graph, partial_paths, current_query_nodes, visited_query_nodes)
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Nodes After Sum: ", visited_query_nodes)
        end

        child_node = popfirst!(node_order)
        parent_idx = 0
        for neighbor in all_neighbors(query_graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
            end
        end
        push!(current_query_nodes, child_node)
        push!(visited_query_nodes, child_node)
        new_partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
        for path in keys(partial_paths)
            running_bounds = partial_paths[path]
            parent_color = only(path[parent_idx])
            for child_color in keys(summary.edgeAvgDegree[parent_color])
                new_path = copy(path)
                push!(new_path, child_color)
                new_bounds = [running_bounds[1]*summary.edgeMinDegree[parent_color][child_color],
                              running_bounds[2]*summary.edgeAvgDegree[parent_color][child_color],
                              running_bounds[3]*summary.edgeMaxDegree[parent_color][child_color],
                ]
                new_partial_paths[new_path] = new_bounds
            end
        end
        partial_paths = new_partial_paths
    end

    final_bounds = [0,0,0]
    for path in keys(partial_paths) 
        final_bounds = final_bounds .+ partial_paths[path]
    end 
    return final_bounds
end

# We use the same general structure to calculate the exact size of the query by finding all paths
# on the original data graph and giving each path a weight of 1.
function get_exact_size(query_graph::DiGraph, data_graph::Graph)
    if is_cyclic(query_graph)
        println("Query Graph Must Be Acyclic")
        return -1
    end

    nodeOrder = topological_sort_by_dfs(query_graph)
    partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
    visited_query_nodes = []
    current_query_nodes = []
    parent_node = popfirst!(nodeOrder)
    child_node = popfirst!(nodeOrder)
    push!(visited_query_nodes, parent_node)
    push!(visited_query_nodes, child_node)
    push!(current_query_nodes, parent_node)
    push!(current_query_nodes, child_node)
    for c1 in vertices(data_graph)
        for c2 in all_neighbors(data_graph, c1)
            partial_paths[[c1,c2]] = [1]
        end
    end

    while length(nodeOrder) > 0
        prev_query_nodes = copy(current_query_nodes)
        for node in prev_query_nodes
            has_living_neighbors = false
            for neighbor in all_neighbors(query_graph, node)
                if !(neighbor in visited_query_nodes)
                    has_living_neighbors = true
                end
            end
            if !has_living_neighbors
                sum_over_node!(partial_paths, current_query_nodes, node)
            end
        end
        
        child_node = popfirst!(nodeOrder)
        parent_idx = 0
        for neighbor in all_neighbors(query_graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
            end
        end
        
        push!(current_query_nodes, child_node)
        push!(visited_query_nodes, child_node)
        new_partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
        for path in keys(partial_paths)
            parent_node = only(path[parent_idx])
            for child_node in all_neighbors(data_graph, parent_node)
                new_path = copy(path)
                push!(new_path, child_node)
                new_partial_paths[new_path] = partial_paths[path]
            end
        end
        partial_paths = new_partial_paths
    end

    final_bounds = [0]
    for path in keys(partial_paths) 
        final_bounds = final_bounds .+ partial_paths[path]
    end 
    return final_bounds
end