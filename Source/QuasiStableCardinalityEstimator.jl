# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
# It currently only handles query graphs without labels.
include("PropertyGraph.jl")
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
using QuasiStableColors

# The ColorSummary struct holds statistical information associated with the colored graph.
# It keeps detailed information about the number of edges between colors of a particular color and which land in
# a particular color. Note that `-1` is used to represent a "wildcard" label. These do not appear in the data graph,
# but they do occur in the query graph.
struct ColorSummary
    color_label_cardinality::Dict{Int, Dict{Int, Int}} # color_label_cardinality[c][v] = num_vertices
    edge_min_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_min_out_deg[c1][c2][e][v2] = min
    edge_min_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_min_in_deg[c1][c2][e][v2] = min
    edge_avg_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_avg_out_deg[c1][c2][e][v2] = avg
    edge_avg_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_avg_in_deg[c1][c2][e][v2] = avg
    edge_max_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_max_out_deg[c1][c2][e][v2] = max
    edge_max_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_max_in_deg[c1][c2][e][v2] = max
end

function generate_color_summary(g::PropertyGraph, numColors::Int)
    color_cardinality = Dict()
    color_label_cardinality = Dict()
    C = q_color(g.graph, n_colors=numColors)
    color_hash::Dict{Int, Int32} = Dict()
    for (color, nodes) in enumerate(C)
        color_label_cardinality[color] = counter(Int)
        color_cardinality[color] = length(nodes)
        for x in nodes
            color_hash[x] = color
            inc!(color_label_cardinality[color], -1)
            for label in g.vertex_labels[x]
                inc!(color_label_cardinality[color], label)
            end
        end
    end

    # We keep separate degree statistics for in-degree and out-degree.
    color_to_color_out_counter::Dict{Int32, Dict{Int32, Any}} = Dict()
    for x in vertices(g.graph)
        c1 = color_hash[x]
        for y in outneighbors(g.graph,x)
            c2 = color_hash[y]
            edge_labels = g.edge_labels[x][y];
            vertex_labels = g.vertex_labels[y];
            if !haskey(color_to_color_out_counter, c1)
                color_to_color_out_counter[c1] = Dict()
            end
            if !haskey(color_to_color_out_counter[c1], c2)
                color_to_color_out_counter[c1][c2] = Dict()
            end
            
            # Since each edge/vertex can have multiple labels associated with it,
            # we must count each edge/vertex label separately in our counter. 
            # Additionally, we need to include a count for the implicit `-1` wildcard label.
            for edge_label in edge_labels
                if !haskey(color_to_color_out_counter[c1][c2], edge_label)
                    color_to_color_out_counter[c1][c2][edge_label] = Dict()
                end
                if !haskey(color_to_color_out_counter[c1][c2][edge_label], -1)
                    color_to_color_out_counter[c1][c2][edge_label][-1] = counter(Int)
                end
                inc!(color_to_color_out_counter[c1][c2][edge_label][-1], x)
                for vertex_label in vertex_labels
                    if !haskey(color_to_color_out_counter[c1][c2][edge_label], vertex_label)
                        color_to_color_out_counter[c1][c2][edge_label][vertex_label] = counter(Int)
                    end
                    inc!(color_to_color_out_counter[c1][c2][edge_label][vertex_label], x)
                end
            end
            if !haskey(color_to_color_out_counter[c1][c2], -1)
                color_to_color_out_counter[c1][c2][-1] = Dict()
            end
            for vertex_label in vertex_labels
                if !haskey(color_to_color_out_counter[c1][c2][-1], vertex_label)
                    color_to_color_out_counter[c1][c2][-1][vertex_label] = counter(Int)
                end
                inc!(color_to_color_out_counter[c1][c2][-1][vertex_label], x)
            end
            if !haskey(color_to_color_out_counter[c1][c2][-1], -1)
                color_to_color_out_counter[c1][c2][-1][-1] = counter(Int)
            end
            inc!(color_to_color_out_counter[c1][c2][-1][-1], x)
        end
    end

    edge_min_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_avg_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_max_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    for c1 in keys(color_to_color_out_counter)
        edge_min_out_deg[c1] = Dict()
        edge_avg_out_deg[c1] = Dict()
        edge_max_out_deg[c1] = Dict()
        for c2 in keys(color_to_color_out_counter[c1])
            edge_min_out_deg[c1][c2] = Dict()
            edge_avg_out_deg[c1][c2] = Dict()
            edge_max_out_deg[c1][c2] = Dict()
            for edge_label in keys(color_to_color_out_counter[c1][c2]) 
                edge_min_out_deg[c1][c2][edge_label] = Dict()
                edge_avg_out_deg[c1][c2][edge_label] = Dict()
                edge_max_out_deg[c1][c2][edge_label] = Dict()
                for vertex_label in keys(color_to_color_out_counter[c1][c2][edge_label])
                    edge_min_out_deg[c1][c2][edge_label][vertex_label] = nv(g.graph) # set this to the max possible value for comparison later
                    edge_max_out_deg[c1][c2][edge_label][vertex_label] = 0 # set to min possible value for comparison later
                    for v in values(color_to_color_out_counter[c1][c2][edge_label][vertex_label])
                        edge_min_out_deg[c1][c2][edge_label][vertex_label] = min(v, edge_min_out_deg[c1][c2][edge_label][vertex_label])
                        edge_max_out_deg[c1][c2][edge_label][vertex_label] = max(v, edge_max_out_deg[c1][c2][edge_label][vertex_label])
                    end
                    edge_avg_out_deg[c1][c2][edge_label][vertex_label] = sum(values(color_to_color_out_counter[c1][c2][edge_label][vertex_label])) / color_cardinality[c1]
                    
                    # if the number of connections is less than the number of vertices in the color,
                    # we can't guarantee the minimum bounds since they won't all map to the same vertex
                    if length(values(color_to_color_out_counter[c1][c2][edge_label][vertex_label])) < color_cardinality[c1]
                        edge_min_out_deg[c1][c2][edge_label][vertex_label] = 0;
                    end
                end
            end
        end
    end

    color_to_color_in_counter::Dict{Int32, Dict{Int32, Any}} = Dict()
    for x in vertices(g.graph)
        c1 = color_hash[x]
        for y in inneighbors(g.graph,x)
            c2 = color_hash[y]
            edge_labels = g.edge_labels[y][x];
            vertex_labels = g.vertex_labels[y];
            if !haskey(color_to_color_in_counter, c1)
                color_to_color_in_counter[c1] = Dict()
            end
            if !haskey(color_to_color_in_counter[c1], c2)
                color_to_color_in_counter[c1][c2] = Dict()
            end
            
            # Since each edge/vertex can have multiple labels associated with it,
            # we must count each edge/vertex label separately in our counter. 
            # Additionally, we need to include a count for the implicit `-1` wildcard label.
            for edge_label in edge_labels
                if !haskey(color_to_color_in_counter[c1][c2], edge_label)
                    color_to_color_in_counter[c1][c2][edge_label] = Dict()
                end
                if !haskey(color_to_color_in_counter[c1][c2][edge_label], -1)
                    color_to_color_in_counter[c1][c2][edge_label][-1] = counter(Int)
                end
                inc!(color_to_color_in_counter[c1][c2][edge_label][-1], x)
                for vertex_label in vertex_labels
                    if !haskey(color_to_color_in_counter[c1][c2][edge_label], vertex_label)
                        color_to_color_in_counter[c1][c2][edge_label][vertex_label] = counter(Int)
                    end
                    inc!(color_to_color_in_counter[c1][c2][edge_label][vertex_label], x)
                end
            end
            if !haskey(color_to_color_in_counter[c1][c2], -1)
                color_to_color_in_counter[c1][c2][-1] = Dict()
            end
            for vertex_label in vertex_labels
                if !haskey(color_to_color_in_counter[c1][c2][-1], vertex_label)
                    color_to_color_in_counter[c1][c2][-1][vertex_label] = counter(Int)
                end
                inc!(color_to_color_in_counter[c1][c2][-1][vertex_label], x)
            end
            if !haskey(color_to_color_in_counter[c1][c2][-1], -1)
                color_to_color_in_counter[c1][c2][-1][-1] = counter(Int)
            end
            inc!(color_to_color_in_counter[c1][c2][-1][-1], x)
        end
    end

    edge_min_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_avg_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_max_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    for c1 in keys(color_to_color_in_counter)
        edge_min_in_deg[c1] = Dict()
        edge_avg_in_deg[c1] = Dict()
        edge_max_in_deg[c1] = Dict()
        for c2 in keys(color_to_color_in_counter[c1])
            edge_min_in_deg[c1][c2] = Dict()
            edge_avg_in_deg[c1][c2] = Dict()
            edge_max_in_deg[c1][c2] = Dict()
            for edge_label in keys(color_to_color_in_counter[c1][c2]) 
                edge_min_in_deg[c1][c2][edge_label] = Dict()
                edge_avg_in_deg[c1][c2][edge_label] = Dict()
                edge_max_in_deg[c1][c2][edge_label] = Dict()
                for vertex_label in keys(color_to_color_in_counter[c1][c2][edge_label])
                    edge_min_in_deg[c1][c2][edge_label][vertex_label] = nv(g.graph) # set this to the max possible value for comparison later
                    edge_max_in_deg[c1][c2][edge_label][vertex_label] = 0 # set to min possible value for comparison later
                    for v in values(color_to_color_in_counter[c1][c2][edge_label][vertex_label])
                        edge_min_in_deg[c1][c2][edge_label][vertex_label] = min(v, edge_min_in_deg[c1][c2][edge_label][vertex_label])
                        edge_max_in_deg[c1][c2][edge_label][vertex_label] = max(v, edge_max_in_deg[c1][c2][edge_label][vertex_label])
                    end
                    edge_avg_in_deg[c1][c2][edge_label][vertex_label] = sum(values(color_to_color_in_counter[c1][c2][edge_label][vertex_label])) / color_cardinality[c1]
                    
                    # if the number of connections is less than the number of vertices in the color,
                    # we can't guarantee the minimum bounds since they won't all map to the same vertex
                    if length(values(color_to_color_in_counter[c1][c2][edge_label][vertex_label])) < color_cardinality[c1]
                        edge_min_in_deg[c1][c2][edge_label][vertex_label] = 0;
                    end
                end
            end
        end
    end
    return ColorSummary(color_label_cardinality, edge_min_out_deg, edge_min_in_deg, 
                                    edge_avg_out_deg, edge_avg_in_deg, edge_max_out_deg, edge_max_in_deg)
end 

function get_color_summary_size(summary)
    numEntries = 0
    for c1 in keys(summary.edge_avg_out_deg)
        for c2 in keys(summary.edge_avg_out_deg[c1])
            for e in keys(summary.edge_avg_out_deg[c1][c2])
                for v in keys(summary.edge_avg_out_deg[c1][c2][e])
                    numEntries += 1
                end
            end
        end
    end
    for c1 in keys(summary.edge_avg_in_deg)
        for c2 in keys(summary.edge_avg_in_deg[c1])
            for e in keys(summary.edge_avg_in_deg[c1][c2])
                for v in keys(summary.edge_avg_in_deg[c1][c2][e])
                    numEntries += 1
                end
            end
        end
    end
    return numEntries*3
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

function sum_over_finished_query_nodes!(query::PropertyGraph, partial_paths, current_query_nodes, visited_query_edges)
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

function get_cardinality_bounds_given_starting_node(query::PropertyGraph, summary::ColorSummary, 
                                                    starting_node::Int; use_partial_sums = true, verbose = false)
    node_order = topological_sort_by_dfs(bfs_tree(Graph(query.graph), starting_node))
    # since we have node labels, should the partial path array have each value instead be subarrays of color-label pairs?
    # Kyle: Because the label is implied by the color -> query_graph_vertex mapping stored in current_query_nodes, I don't think
    # we have to keep the label in the partial paths object.
    partial_paths::Array{Tuple{Array{Int}, Array{Float64}}} = [] # each tuple contains a pairing of color paths -> bounds
    visited_query_edges = []
    current_query_nodes = []
    # Change this to initialize the partial_paths as 1-node paths with bounds corresponding to the 
    # number of nodes of the particular label within each color

    # new initialization 
    parent_node = popfirst!(node_order)
    parent_label = query.vertex_labels[parent_node][1]
    push!(current_query_nodes, parent_node)
    # initialize partial_paths with all possible starting color/vertex possibilities.
    for color in keys(summary.color_label_cardinality)
        # only use the parent label
        if (haskey(summary.color_label_cardinality[color], parent_label))
            push!(partial_paths, ([color], [summary.color_label_cardinality[color][parent_label],
                                                            summary.color_label_cardinality[color][parent_label],
                                                            summary.color_label_cardinality[color][parent_label]]))
        end
    end
    child_node = parent_node
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

        # get the next child from the node order
        child_node = popfirst!(node_order)
        parent_idx = 0
        outEdge = false
        for neighbor in all_neighbors(query.graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                # gets the index of the parent in the list of current nodes
                parent_idx = indexin(neighbor, current_query_nodes)
                if parent_node in inneighbors(query.graph, child_node)
                    outEdge = true
                end
                break
            end
        end
        # push the current edge and nodes to the visited lists
        if outEdge 
            push!(visited_query_edges, (parent_node, child_node))
        else 
            push!(visited_query_edges, (child_node, parent_node))
        end
        push!(current_query_nodes, child_node)

        # get the appropriate labels, the query only uses one label per vertex/node
        edge_label = 0
        if outEdge
            edge_label = query.edge_labels[parent_node][child_node][1]
        else
            edge_label =  query.edge_labels[child_node][parent_node][1]
        end 
        child_label = query.vertex_labels[child_node][1]

        # update the partial paths using the parent-child combo that comes next from the query
        new_partial_paths::Array{Tuple{Array{Int}, Array{Float64}}} = []
        for path_and_bounds in partial_paths
            path = path_and_bounds[1] # using a tuple causes intermediate data structure?
            running_bounds = path_and_bounds[2]
            parent_color = only(path[parent_idx])
            # account for colors with no outgoing children
            if outEdge && haskey(summary.edge_avg_out_deg, parent_color)
                for child_color in keys(summary.edge_avg_out_deg[parent_color])
                    new_path = copy(path)
                    push!(new_path, child_color)     
                    if (haskey(summary.edge_avg_out_deg[parent_color][child_color], edge_label)) &&
                            (haskey(summary.edge_avg_out_deg[parent_color][child_color][edge_label], child_label))
                        new_bounds = [running_bounds[1]*summary.edge_min_out_deg[parent_color][child_color][edge_label][child_label],
                            running_bounds[2]*summary.edge_avg_out_deg[parent_color][child_color][edge_label][child_label],
                            running_bounds[3]*summary.edge_max_out_deg[parent_color][child_color][edge_label][child_label],
                        ]
                        push!(new_partial_paths, (new_path, new_bounds))
                    end 
                end
            elseif !outEdge && haskey(summary.edge_avg_in_deg, parent_color)
                for child_color in keys(summary.edge_avg_in_deg[parent_color])
                    new_path = copy(path)
                    push!(new_path, child_color)     
                    if (haskey(summary.edge_avg_in_deg[parent_color][child_color], edge_label)) &&
                            (haskey(summary.edge_avg_in_deg[parent_color][child_color][edge_label], child_label))
                        new_bounds = [running_bounds[1]*summary.edge_min_in_deg[parent_color][child_color][edge_label][child_label],
                            running_bounds[2]*summary.edge_avg_in_deg[parent_color][child_color][edge_label][child_label],
                            running_bounds[3]*summary.edge_max_in_deg[parent_color][child_color][edge_label][child_label],
                        ]
                        push!(new_partial_paths, (new_path, new_bounds))
                    end
                end
            end
        end
        partial_paths = new_partial_paths
    end

    # To account for cyclic queries, we check whether there are any remaining edges that have not
    # been processed. If so, we set the lower bound to 0, reduce the average estimate accordingly, and leave
    # the upper bound unchanged.
    remaining_edges = []
    for edge in edges(query.graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end

    # Sum over the calculated partial paths to get the final bounds
    final_bounds = [0,0,0]
    for path_and_bounds in partial_paths
        path = path_and_bounds[1]
        bounds = path_and_bounds[2]
        lower = only(bounds[1])
        if length(remaining_edges) > 0
            lower = 0
        end
        average = only(bounds[2])
        # scale down the average if there are remaining non-tree-edges
        for edge in remaining_edges
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_color = only(path[parent_node_idx])
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_color = only(path[child_node_idx][1])
            child_label = only(path[child_node_idx][1])
            edge_label = query.edge_labels[edge[1]][edge[2]][1]
            probability_of_edge = 0
            if (haskey(summary.edge_avg_out_deg, parent_color) 
                    && haskey(summary.edge_avg_out_deg[parent_color], child_color)
                        && haskey(summary.edge_avg_out_deg[parent_color][child_color], edge_label)
                            && haskey(summary.edge_avg_out_deg[parent_color][child_color][edge_label], child_label))
                probability_of_edge = summary.edge_avg_out_deg[parent_color][child_color][edge_label][child_label]/summary.color_label_cardinality[child_color][child_label]
            end
            average *= probability_of_edge
        end
        upper = only(bounds[3])
        final_bounds = final_bounds .+ [lower, average, upper]
    end
    return final_bounds
end

function get_cardinality_bounds(query::PropertyGraph, summary::ColorSummary; 
                                use_partial_sums = true, try_all_starting_nodes=true, verbose = false)
    if try_all_starting_nodes
        final_bounds = [0, 0, Inf]
        for node in vertices(query.graph)
            cur_bounds = get_cardinality_bounds_given_starting_node(query, summary, node;
                                                                    use_partial_sums=use_partial_sums, verbose=verbose)
            final_bounds[1] = max(final_bounds[1], cur_bounds[1])
            final_bounds[2] += cur_bounds[2]
            final_bounds[3] = min(final_bounds[3], cur_bounds[3])
        end
        final_bounds[2] /= nv(query.graph)
        return final_bounds
    end
    return get_cardinality_bounds_given_starting_node(query, summary, vertices(query.graph)[1];
                                                      use_partial_sums=use_partial_sums, verbose=verbose)
end

# We use the same general structure to calculate the exact size of the query by finding all paths
# on the original data graph and giving each path a weight of 1. 
function get_exact_size(query::PropertyGraph, data::PropertyGraph; use_partial_sums = true, verbose=false)
    node_order = topological_sort_by_dfs(bfs_tree(Graph(query.graph), vertices(query.graph)[1]))
    # right now the path is an array of node-label tuples, but the "label" part isn't necessary since
    # we have access to the original data graph. It's there because the summing function needs it?
    partial_paths::Array{Tuple{Array{Int}, Int}} = []
    visited_query_edges = []
    current_query_nodes = []
    if verbose
        println("Node Order: ", node_order)
    end
    # initialize partial_paths as 1-node paths with label matching the
    # initial query node's label
    parent_node = popfirst!(node_order)
    child_node = parent_node
    parent_label = query.vertex_labels[parent_node][1]
    push!(current_query_nodes, parent_node)
    for node in vertices(data.graph)
        node_labels = data.vertex_labels[node]
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
        if use_partial_sums
            sum_over_finished_query_nodes!(query, partial_paths, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end

        child_node = popfirst!(node_order)
        parent_idx = 0
        outEdge = false
        for neighbor in all_neighbors(query.graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                parent_idx = indexin(neighbor, current_query_nodes)
                if neighbor in inneighbors(query.graph, child_node)
                    outEdge = true
                end
                break
            end
        end
        query_edge_label = 0
        if outEdge
            query_edge_label = query.edge_labels[parent_node][child_node][1]
        else
            query_edge_label =  query.edge_labels[child_node][parent_node][1]
        end 
        query_child_label = query.vertex_labels[child_node][1]
        
        push!(current_query_nodes, child_node)
        push!(visited_query_edges, (parent_node, child_node))
        new_partial_paths::Array{Tuple{Array{Int}, Int}} = []
        for path_and_weight in partial_paths
            path = path_and_weight[1]
            weight = path_and_weight[2]
            parent_node = only(path[parent_idx])
            potential_child_nodes = []
            if outEdge
                for data_child_node in outneighbors(data.graph, parent_node)
                    new_weight = weight
                    # only add a new partial path if the edge label and node label match our query
                    # do we not have to worry about matching the parent node?
                    data_edge_labels = data.edge_labels[parent_node][data_child_node]
                    data_child_labels = data.vertex_labels[data_child_node]
                    if (query_child_label == -1  || in(query_child_label, data_child_labels)) && 
                        (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                        new_path = copy(path)
                        push!(new_path, data_child_node)
                        push!(new_partial_paths, (new_path, new_weight))
                    end
                end
            else
                for data_child_node in inneighbors(data.graph, parent_node)
                    new_weight = weight
                    # only add a new partial path if the edge label and node label match our query
                    # do we not have to worry about matching the parent node?
                    data_edge_labels = data.edge_labels[data_child_node][parent_node]
                    data_child_labels = data.vertex_labels[data_child_node]
                    if (query_child_label == -1  || in(query_child_label, data_child_labels)) && 
                        (query_edge_label == -1 || in(query_edge_label, data_edge_labels))
                        new_path = copy(path)
                        push!(new_path, data_child_node)
                        push!(new_partial_paths, (new_path, new_weight))
                    end
                end
            end
        end
        partial_paths = new_partial_paths
    end
    remaining_edges = []
    for edge in edges(query.graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end
    exact_size = 0
    for path_and_weight in partial_paths 
        path = path_and_weight[1]
        weight = path_and_weight[2]
        satisfies_cycles = true
        for edge in remaining_edges
            # Only count the cycle as satisfied if this remaining edge's label matches the query graph's edge label
            # how do we check this if the remaining edges are taken from the query graph? Isn't this guaranteed?

            # get the parent node from the list of current query nodes
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_data_node = only(path[parent_node_idx])
            # get the child node from the list of current query nodes
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_data_node = only(path[child_node_idx])

            # check if the edge label exists, if it doesn't then we can break here
            # don't need to check parent node because we got the parent node from the data graph,
            # but we do need to check if there is an edge connection to the child
            if (!haskey(data.edge_labels[parent_data_node], child_data_node))
                satisfies_cycles = false;
                break;
            end
            data_edge_labels = data.edge_labels[parent_data_node][child_data_node]
            data_child_vertex_labels = data.vertex_labels[child_data_node]
            query_edge_label = only(query.edge_labels[edge[1]][edge[2]])
            query_child_vertex_label = only(query.vertex_labels[edge[2]])
            if !((query_edge_label == -1 || in(query_edge_label, data_edge_labels)) && 
                    (query_child_vertex_label == -1 || in(query_child_vertex_label, data_child_vertex_labels)))
                satisfies_cycles = false
                break
            end
            if query_edge_label == -1
                weight *= length(data_edge_labels)
            end
        end 
        if satisfies_cycles
            exact_size = exact_size + weight
        end
    end 
    return exact_size
end
