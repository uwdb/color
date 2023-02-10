
# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
# It currently only handles query graphs without labels.
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
using QuasiStableColors

# Create a property graph struct which has the following things: underlying graph,
# and dictionaries which map edges to label sets and vertices to label sets
struct PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Int, Dict{Int, Set{Int}}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Dict{Int, Set{Int}} # vertex_labels[n] = { labels }
end

# Everywhere that we reference a graph (query or data) we should replace it with our new struct.


# Extend the colorsummary to account for edge labels and vertex labels as additional layers of dict nesting
# Change the color_cardinality to account for node labels and remove the edge_cardinality dict
struct ColorSummary
    color_cardinality::Dict{Int, Dict{Int, Int}} # color_cardinality[c][v] = num_vertices
    edge_min_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_min_deg[c1][c2][e][v2] = min
    edge_avg_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_avg_deg[c1][c2][v2] = avg
    edge_max_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_max_deg[c1][c2][v2] = max
end

function generate_color_summary(g::PropertyGraph, numColors::Int)
    color_cardinality = Dict()
    undirected_graph = Graph(g.graph)
    C = q_color(undirected_graph, n_colors=numColors)
    color_hash::Dict{Int, Int} = Dict()
    for (color, nodes) in enumerate(C)
        for x in nodes
            color_hash[x] = color
            color_cardinality[color] = counter()
            # for each kind of label, keep track 
            for label in g.vertex_labels[x]
                inc!(color_cardinality[label], label)
            end
        end
    end

    color_to_color_counter::Dict{Int, Dict{Int, Any}} = Dict()
    for x in vertices(g)
        c1 = color_hash[x]
        for y in outneighbors(g.graph,x) 
            c2 = color_hash[y]
            edge_label = g.edge_labels[x][y];
            vertex_label = g.vertex_labels[y];
            if !haskey(color_to_color_counter, c1)
                color_to_color_counter[c1] = Dict()
            end
            if !haskey(color_to_color_counter[c1], c2)
                color_to_color_counter[c1][c2] = Dict()
            end
            if !haskey(color_to_color_counter[c1][c2], edge_label)
                color_to_color_counter[c1][c2][edge_label] = Dict()
            end
            if !haskey(color_to_color_counter[c1][c2][edge_label], vertex_label)
                color_to_color_counter[c1][c2][edge_label][vertex_label] = counter(Int)
            end
            inc!(color_to_color_counter[c1][c2][edge_label][vertex_label], x)
        end
    end

    # edge_cardinality::Dict{Int, Dict{Int, Float64}} = Dict()
    edge_min_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_avg_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_max_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    for c1 in keys(color_to_color_counter)
        edge_min_deg[c1] = Dict()
        edge_avg_deg[c1] = Dict()
        edge_max_deg[c1] = Dict()
        for c2 in keys(color_to_color_counter[c1])
            edge_min_deg[c1][c2] = Dict()
            edge_avg_deg[c1][c2] = Dict()
            edge_max_deg[c1][c2] = Dict()
            for edge_label in keys(color_to_color_counter[c1][c2]) 
                edge_min_deg[c1][c2][edge_label] = Dict()
                edge_avg_deg[c1][c2][edge_label] = Dict()
                edge_max_deg[c1][c2][edge_label] = Dict()
                # not sure what to do for averages yet
                for vertex_label in keys(color_to_color_counter[c1][c2][edge_label])
                    edge_min_deg[c1][c2][edge_label][vertex_label] = nv(g.graph) # set this to the max possible value for comparison later
                    edge_max_deg[c1][c2][edge_label][vertex_label] = 0 # set to min possible value for comparison later
                    for v in values(color_to_color_counter[c1][c2][edge_label][vertex_label])
                        edge_min_deg[c1][c2][edge_label][vertex_label] = min(v, edge_min_deg[c1][c2][edge_label][vertex_label])
                        edge_max_deg[c1][c2][edge_label][vertex_label] = max(v, edge_max_deg[c1][c2][edge_label][vertex_label])
                    end

                    # previously, to get edge_avg_deg divided the edge_cardinality by color cardinality
                    # edge_cardinality represents the number of edges from c1 to c2
                    # color_cardinality represents the number of nodes in c1
                    # so now would we use the color_to_color_counter instead of edge_cardinality?
                    edge_avg_deg[c1][c2][edge_label][vertex_label] = sum(values(color_to_color_counter[c1][c2][edge_label][vertex_label])) / color_cardinality[c1];
                    # if the number of connections is less than the number of vertices in the color,
                    # we can't guarantee the minimum bounds since they won't all map to the same vertex
                    if length(values(color_to_color_counter[c1][c2][edge_label][vertex_label])) < color_cardinality[c1]
                        edge_min_deg[c1][c2][edge_label][vertex_label] = 0;
                    end
                end
            end
        end
    end

    return ColorSummary(color_cardinality, edge_min_deg, edge_avg_deg, edge_max_deg)
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

function get_cardinality_bounds_given_starting_node(query::PropertyGraph, summary::ColorSummary, 
                                                    starting_node::Int; use_partial_sums = true, verbose = false)
    node_order = topological_sort_by_dfs(bfs_tree(query.graph, starting_node))
    # since we have node labels, should the partial path array have each value instead be subarrays of color-label pairs?
    partial_paths::Dict{Array{Array{Int}}, Array{Float64}} = Dict()
    visited_query_edges = []
    current_query_nodes = []
    # Change this to initalize the partial_paths as 1-node paths with bounds corresponding to the 
    # number of nodes of the particular label within each color

    # new initialization:
    parent_node = popfirst!(node_order)
    push!(current_query_nodes, parent_node)
    for color in keys(summary.color_cardinality)
        for vertex_label in keys(summary.color_cardinality[color])
            partial_paths[[[color, vertex_label]]] = [summary.color_cardinality[color][vertex_label],
                                                        summary.color_cardinality[color][vertex_label],
                                                        summary.color_cardinality[color][vertex_label]]
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
        for neighbor in inneighbors(query.graph, child_node)
            if neighbor in current_query_nodes
                parent_node = neighbor
                # gets the index of the parent in the list of current nodes
                parent_idx = indexin(neighbor, current_query_nodes)
            end
        end
        # push the current edge and nodes to the visited lists
        # if we want information about the labels we can refer to the property graph
        push!(visited_query_edges, (parent_node, child_node))
        push!(current_query_nodes, child_node)

        # update the partial paths using the parent-child combo that comes next from the query
        new_partial_paths::Dict{Array{Array{Int}}, Array{Float64}} = Dict()
        for path in keys(partial_paths)
            running_bounds = partial_paths[path]
            parent_info = only(path[parent_idx])
            parent_color = parent_info[1]
            for child_color in keys(summary.edge_avg_deg[parent_color])
                # are these the right labels to use?
                for edge_label in keys(summary.edge_avg_deg[parent_color][child_color])
                    for child_vertex_label in keys(summary.edge_avg_deg[parent_color][child_color][edge_label])
                        new_path = copy(path)
                        push!(new_path, (child_color, child_vertex_label))
                        new_bounds = [running_bounds[1] * summary.edge_min_deg[parent_color][child_color][edge_label][child_vertex_label],
                                        running_bounds[1] * summary.edge_avg_deg[parent_color][child_color][edge_label][child_vertex_label],
                                        running_bounds[1] * summary.edge_max_deg[parent_color][child_color][edge_label][child_vertex_label]]
                        new_partial_paths[new_path] = new_bounds
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
    for edge in edges(query_graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end

    final_bounds = [0,0,0]
    for path in keys(partial_paths) 
        lower = only(partial_paths[path][1])
        if length(remaining_edges) > 0
            lower = 0
        end
        average = only(partial_paths[path][2])
        for edge in remaining_edges
            probability_of_edge = 0
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_info = only(path[parent_node_idx])
            parent_color = parent_info[1]
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_info = only(path[child_node_idx])
            child_color = child_info[1]
            child_label = child_info[2]
            # get all of the labels for the current edge
            for edge_label in query.edge_labels[edge[1]][edge[2]]
                if haskey(summary.edge_avg_deg, parent_color) & haskey(summary.edge_avg_deg[parent_color], child_color)
                & haskey(summary.edge_avg_deg[parent_color][child_color], edge_label)
                & haskey(summary.edge_avg_deg[parent_color][child_color][edge_label], child_label)
                    probability_of_edge = summary.edge_avg_deg[parent_color][child_color][edge_label][child_label]/summary.color_cardinality[child_color]
                end
                average *= probability_of_edge
            end
        end
        
        upper = only(partial_paths[path][3])
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

#Make similar changes to get exact size
function get_exact_size(query_graph::DiGraph, data_graph::DiGraph; use_partial_sums = true, verbose=false)
    node_order = topological_sort_by_dfs(bfs_tree(query_graph, vertices(query_graph)[1]))
    partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
    visited_query_edges = []
    current_query_nodes = []
    parent_node = popfirst!(node_order)
    child_node = popfirst!(node_order)
    push!(visited_query_edges, (parent_node, child_node))
    push!(current_query_nodes, parent_node)
    push!(current_query_nodes, child_node)
    for c1 in vertices(data_graph)
        for c2 in outneighbors(data_graph, c1)
            partial_paths[[c1,c2]] = [1]
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
        new_partial_paths::Dict{Array{Int}, Array{Float64}} = Dict()
        for path in keys(partial_paths)
            parent_node = only(path[parent_idx])
            for child_node in outneighbors(data_graph, parent_node)
                # only add a new partial path if the edge label and node label match our query
                new_path = copy(path)
                push!(new_path, child_node)
                new_partial_paths[new_path] = partial_paths[path]
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

    final_bounds = [0]
    for path in keys(partial_paths) 
        satisfies_cycles = true
        for edge in remaining_edges
            # Only count the cycle as satisfied if this remaining edge's label matches the query graph's edge label
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_data_node = only(path[parent_node_idx])
            child_node_idx = indexin(edge[2], current_query_nodes)
            child_data_node = only(path[child_node_idx])
            if !(child_data_node in inneighbors(data_graph, parent_data_node))
                satisfies_cycles = false
            end
        end 
        if satisfies_cycles
            final_bounds = final_bounds .+ partial_paths[path]
        end
    end 
    return final_bounds
end
