# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.
# It currently only handles query graphs without labels.
include("PropertyGraph.jl")
using Graphs
using QuasiStableColors
using Probably

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
    color_filters::Dict{Int, BloomFilter} # color_filters[c] = filter
end

function generate_color_summary(g::DataGraph, numColors::Int; weighting=true, verbose=false)
    QSC = QuasiStableColors
    color_filters = Dict()
    color_cardinality = Dict()
    color_label_cardinality = Dict()
    if (verbose) 
        println("Started coloring")
    end
    C = QSC.q_color(g.graph, n_colors=numColors, weighting=weighting)
    if (verbose)
        println("Finished coloring")
    end
    color_hash = QSC.node_map(C)

    # initialize color filters for data labels
    current_color = 1;
    if (verbose)
        println("Started bloom filters")
    end
    for color in C.partitions
        num_nodes = only(size(color))
        accepted_error = 0.01
        parameters = constrain(BloomFilter, fpr=accepted_error, capacity=num_nodes)
        color_filters[current_color] = BloomFilter(parameters.m, parameters.k)
        current_color += 1
    end
    if (verbose)
        println("Finished bloom filters")
    end

    if (verbose)
        println("Started cardinality counts")
    end
    for node in keys(color_hash)
        color = color_hash[node]
        data_label = get_data_label(g, node)
        if (data_label != -1)
            push!(color_filters[color], data_label)
        end

        # initialize color-label cardinality counter
        if (!haskey(color_label_cardinality, color))
            color_label_cardinality[color] = counter(Int)
        end
        # initialize color cardinality counter
        if (!haskey(color_cardinality, color))
            color_cardinality[color] = 1
        else
            color_cardinality[color] += 1
        end

        # increment counter for all labels, including wildcards
        inc!(color_label_cardinality[color], -1)
        for label in g.vertex_labels[node]
            inc!(color_label_cardinality[color], label)
        end
    end
    if (verbose)
        println("Finished cardinality counts")
    end
    
    if (verbose)
        println("Started tracking statistics")
    end
    # We keep separate degree statistics for in-degree and out-degree.
    color_to_color_out_counter::Dict{Int32, Dict{Int32, Any}} = Dict()
    for x in vertices(g.graph)
        c1 = color_hash[x]
        for y in outneighbors(g.graph,x)
            c2 = color_hash[y]
            edge_labels = []
            copy!(edge_labels, g.edge_labels[(x, y)])
            push!(edge_labels, -1)
            vertex_labels = []
            copy!(vertex_labels, g.vertex_labels[y])
            push!(vertex_labels, -1)
            
            # Since each edge/vertex can have multiple labels associated with it,
            # we must count each edge/vertex label separately in our counter. 
            # Additionally, we need to include a count for the implicit `-1` wildcard label.
            for edge_label in edge_labels
                if !haskey(color_to_color_out_counter, edge_label)
                    color_to_color_out_counter[edge_label] = Dict()
                end
                for vertex_label in vertex_labels
                    if !haskey(color_to_color_out_counter[edge_label], vertex_label)
                        color_to_color_out_counter[edge_label][vertex_label] = Dict()
                    end
                    if !haskey(color_to_color_out_counter[edge_label][vertex_label], c1)
                        color_to_color_out_counter[edge_label][vertex_label][c1] = Dict()
                    end
                    if !haskey(color_to_color_out_counter[edge_label][vertex_label][c1], c2)
                        color_to_color_out_counter[edge_label][vertex_label][c1][c2] = counter(Int)
                    end
                    inc!(color_to_color_out_counter[edge_label][vertex_label][c1][c2], x)
                end
            end
        end
    end

    edge_min_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_avg_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_max_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    for edge_label in keys(color_to_color_out_counter)
        edge_min_out_deg[edge_label] = Dict()
        edge_avg_out_deg[edge_label] = Dict()
        edge_max_out_deg[edge_label] = Dict()
        for vertex_label in keys(color_to_color_out_counter[edge_label])
            edge_min_out_deg[edge_label][vertex_label] = Dict()
            edge_avg_out_deg[edge_label][vertex_label] = Dict()
            edge_max_out_deg[edge_label][vertex_label] = Dict()
            for c1 in keys(color_to_color_out_counter[edge_label][vertex_label]) 
                edge_min_out_deg[edge_label][vertex_label][c1] = Dict()
                edge_avg_out_deg[edge_label][vertex_label][c1] = Dict()
                edge_max_out_deg[edge_label][vertex_label][c1] = Dict()
                for c2 in keys(color_to_color_out_counter[edge_label][vertex_label][c1])
                    edge_min_out_deg[edge_label][vertex_label][c1][c2] = nv(g.graph) # set this to the max possible value for comparison later
                    edge_max_out_deg[edge_label][vertex_label][c1][c2] = 0 # set to min possible value for comparison later
                    for v in values(color_to_color_out_counter[edge_label][vertex_label][c1][c2])
                        edge_min_out_deg[edge_label][vertex_label][c1][c2] = min(v, edge_min_out_deg[edge_label][vertex_label][c1][c2])
                        edge_max_out_deg[edge_label][vertex_label][c1][c2] = max(v, edge_max_out_deg[edge_label][vertex_label][c1][c2])
                    end
                    edge_avg_out_deg[edge_label][vertex_label][c1][c2] = sum(values(color_to_color_out_counter[edge_label][vertex_label][c1][c2])) / color_cardinality[c1]
                    
                    # if the number of connections is less than the number of vertices in the color,
                    # we can't guarantee the minimum bounds since they won't all map to the same vertex
                    if length(values(color_to_color_out_counter[edge_label][vertex_label][c1][c2])) < color_cardinality[c1]
                        edge_min_out_deg[edge_label][vertex_label][c1][c2] = 0;
                    end
                end
            end
        end
    end

    # We keep separate degree statistics for in-degree and out-degree.
    color_to_color_in_counter::Dict{Int32, Dict{Int32, Any}} = Dict()
    for x in vertices(g.graph)
        c1 = color_hash[x]
        for y in inneighbors(g.graph,x)
            c2 = color_hash[y]
            edge_labels = []
            copy!(edge_labels, g.edge_labels[(y, x)])
            push!(edge_labels, -1)
            vertex_labels = []
            copy!(vertex_labels, g.vertex_labels[y])
            push!(vertex_labels, -1)
            
            # Since each edge/vertex can have multiple labels associated with it,
            # we must count each edge/vertex label separately in our counter. 
            # Additionally, we need to include a count for the implicit `-1` wildcard label.
            for edge_label in edge_labels
                if !haskey(color_to_color_in_counter, edge_label)
                    color_to_color_in_counter[edge_label] = Dict()
                end
                for vertex_label in vertex_labels
                    if !haskey(color_to_color_in_counter[edge_label], vertex_label)
                        color_to_color_in_counter[edge_label][vertex_label] = Dict()
                    end
                    if !haskey(color_to_color_in_counter[edge_label][vertex_label], c1)
                        color_to_color_in_counter[edge_label][vertex_label][c1] = Dict()
                    end
                    if !haskey(color_to_color_in_counter[edge_label][vertex_label][c1], c2)
                        color_to_color_in_counter[edge_label][vertex_label][c1][c2] = counter(Int)
                    end
                    inc!(color_to_color_in_counter[edge_label][vertex_label][c1][c2], x)
                end
            end
        end
    end

    edge_min_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_avg_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    edge_max_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} = Dict()
    for edge_label in keys(color_to_color_in_counter)
        edge_min_in_deg[edge_label] = Dict()
        edge_avg_in_deg[edge_label] = Dict()
        edge_max_in_deg[edge_label] = Dict()
        for vertex_label in keys(color_to_color_in_counter[edge_label])
            edge_min_in_deg[edge_label][vertex_label] = Dict()
            edge_avg_in_deg[edge_label][vertex_label] = Dict()
            edge_max_in_deg[edge_label][vertex_label] = Dict()
            for c1 in keys(color_to_color_in_counter[edge_label][vertex_label]) 
                edge_min_in_deg[edge_label][vertex_label][c1] = Dict()
                edge_avg_in_deg[edge_label][vertex_label][c1] = Dict()
                edge_max_in_deg[edge_label][vertex_label][c1] = Dict()
                for c2 in keys(color_to_color_in_counter[edge_label][vertex_label][c1])
                    edge_min_in_deg[edge_label][vertex_label][c1][c2] = nv(g.graph) # set this to the max possible value for comparison later
                    edge_max_in_deg[edge_label][vertex_label][c1][c2] = 0 # set to min possible value for comparison later
                    for v in values(color_to_color_in_counter[edge_label][vertex_label][c1][c2])
                        edge_min_in_deg[edge_label][vertex_label][c1][c2] = min(v, edge_min_in_deg[edge_label][vertex_label][c1][c2])
                        edge_max_in_deg[edge_label][vertex_label][c1][c2] = max(v, edge_max_in_deg[edge_label][vertex_label][c1][c2])
                    end
                    edge_avg_in_deg[edge_label][vertex_label][c1][c2] = sum(values(color_to_color_in_counter[edge_label][vertex_label][c1][c2])) / color_cardinality[c1]
                    
                    # if the number of connections is less than the number of vertices in the color,
                    # we can't guarantee the minimum bounds since they won't all map to the same vertex
                    if length(values(color_to_color_in_counter[edge_label][vertex_label][c1][c2])) < color_cardinality[c1]
                        edge_min_in_deg[edge_label][vertex_label][c1][c2] = 0;
                    end
                end
            end
        end
    end
    if (verbose)
        println("Finished tracking statistics")
    end
    color_to_color_in_counter = Dict()
    color_to_color_out_counter = Dict()
    return ColorSummary(color_label_cardinality, edge_min_out_deg, edge_min_in_deg, 
                                    edge_avg_out_deg, edge_avg_in_deg, edge_max_out_deg, edge_max_in_deg, color_filters)
end

function get_color_summary_size(summary)
    numEntries = 0
    for e in keys(summary.edge_avg_out_deg)
        for v in keys(summary.edge_avg_out_deg[e])
            for c1 in keys(summary.edge_avg_out_deg[e][v])
                for c2 in keys(summary.edge_avg_out_deg[e][v][c1])
                    numEntries += 1
                end
            end
        end
    end
    for e in keys(summary.edge_avg_in_deg)
        for v in keys(summary.edge_avg_in_deg[e])
            for c1 in keys(summary.edge_avg_in_deg[e][v])
                for c2 in keys(summary.edge_avg_in_deg[e][v][c1])
                    numEntries += 1
                end
            end
        end
    end
    return numEntries*3 # To account for min, avg, and max stats
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

function handle_extra_edges!(query::QueryGraph, summary::ColorSummary, partial_paths, current_query_nodes, visited_query_edges)
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

    # Sum over the calculated partial paths to get the final bounds
    for i  in range(1, length(partial_paths))
        path = partial_paths[i][1]
        bounds = partial_paths[i][2]
        lower = only(bounds[1])
        if length(remaining_edges) > 0
            lower = 0
        end
        average = only(bounds[2])
        # scale down the average if there are remaining non-tree-edges
        for edge in remaining_edges
            parent_node_idx = indexin(edge[1], current_query_nodes)
            parent_color = only(path[parent_node_idx])
            new_node_idx = indexin(edge[2], current_query_nodes)
            child_color = only(path[new_node_idx][1])
            child_label = only(query.vertex_labels[edge[2]])
            # don't have to check data label because these nodes are already in the
            # partial path, so we have already ensured that the colors are appropriate
            edge_label = only(query.edge_labels[(edge[1],edge[2])])
            probability_of_edge = 0
            if (haskey(summary.edge_avg_out_deg, edge_label) 
                    && haskey(summary.edge_avg_out_deg[edge_label], child_label) # so we know that the child label is not appearing in the edge label table...
                        && haskey(summary.edge_avg_out_deg[edge_label][child_label], parent_color)
                            && haskey(summary.edge_avg_out_deg[edge_label][child_label][parent_color], child_color))
                probability_of_edge = summary.edge_avg_out_deg[edge_label][child_label][parent_color][child_color]/summary.color_label_cardinality[child_color][child_label]
            end
            average *= probability_of_edge
        end
        upper = only(bounds[3])
        partial_paths[i] = (path, [lower, average, upper])
    end 
end

function sum_over_finished_query_nodes!(query::QueryGraph, partial_paths, current_query_nodes, visited_query_edges)
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

function get_cardinality_bounds(query::QueryGraph, summary::ColorSummary; use_partial_sums = true, verbose = false)
    node_order = get_min_width_node_order(query.graph)
    if verbose
        println("Node Order:", node_order)
    end
    # Because the label is implied by the color -> query_graph_vertex mapping stored in current_query_nodes,
    # we don't have to keep the label in the partial paths object.
    partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}} = [] # each tuple contains a pairing of color paths -> bounds
    visited_query_edges = []
    current_query_nodes = []
    # Change this to initialize the partial_paths as 1-node paths with bounds corresponding to the 
    # number of nodes of the particular label within each color.

    old_node = popfirst!(node_order)
    parent_label = query.vertex_labels[old_node][1]
    parent_data_label = get_data_label(query, old_node)
    push!(current_query_nodes, old_node)
    # Initialize partial_paths with all possible starting color/vertex possibilities.
    for color in keys(summary.color_label_cardinality)
        # Only use the parent label.
        if (haskey(summary.color_label_cardinality[color], parent_label))
            # if the parent has a specified data label, only use colors that the filters approve
            if (parent_data_label == -1 || parent_data_label in summary.color_filters[color])
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
                # Gets the index of the parent in the list of current nodes.
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
        new_data_label = get_data_label(query, new_node)

        # Update the partial paths using the parent-child combo that comes next from the query.
        new_partial_paths::Vector{Tuple{Vector{Int}, Vector{Float64}}} = []
        for path_and_bounds in partial_paths
            path = path_and_bounds[1]
            running_bounds = path_and_bounds[2]
            old_color = only(path[parent_idx])
            # Account for colors with no outgoing children.
            if outEdge && haskey(summary.edge_avg_out_deg, edge_label) && 
                haskey(summary.edge_avg_out_deg[edge_label], new_label) && 
                haskey(summary.edge_avg_out_deg[edge_label][new_label], old_color)
                for new_color in keys(summary.edge_avg_out_deg[edge_label][new_label][old_color])
                    if (new_data_label != -1 && !(new_data_label in summary.color_filters[new_color]))
                        # if the child has a data label that isn't expected to be in this color,
                        # then skip adding this color to the partial paths
                        continue
                    end
                    new_path = copy(path)
                    push!(new_path, new_color)   
                    new_bounds = [running_bounds[1]*summary.edge_min_out_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[2]*summary.edge_avg_out_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[3]*summary.edge_max_out_deg[edge_label][new_label][old_color][new_color],
                                    ]
                    if (new_data_label != -1)
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
                    if (new_data_label != -1 && !(new_data_label in summary.color_filters[new_color]))
                        # if the child has a data label that isn't expected to be in this color,
                        # then skip adding this color to the partial paths
                        continue
                    end                
                    new_path = copy(path)
                    push!(new_path, new_color)
                    new_bounds = [running_bounds[1]*summary.edge_min_in_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[2]*summary.edge_avg_in_deg[edge_label][new_label][old_color][new_color],
                                    running_bounds[3]*summary.edge_max_in_deg[edge_label][new_label][old_color][new_color],
                                    ]
                    if (new_data_label != -1)
                        # we have already confirmed that the data label is in the color, but if the data label isn't -1
                        # then we need to scale down the result since we only want to consider one of the many nodes in the new color
                        new_bounds[2] = new_bounds[2] / summary.color_label_cardinality[new_color][new_label] # is this wrong because it divides by only a part of the color cardinality instead of the whole? still shouldn't result in 0 tho
                        # we also need to set the minimum to 0 but keep the maximum the same
                        new_bounds[1] = 0
                    end
                    push!(new_partial_paths, (new_path, new_bounds))
                end
            end
        end
        partial_paths = new_partial_paths
        handle_extra_edges!(query, summary, partial_paths, current_query_nodes, visited_query_edges)
    end

    # Sum over the calculated partial paths to get the final bounds.
    final_bounds = [0,0,0]
    numZeroAvg = 0
    for path_and_bounds in partial_paths
        final_bounds = final_bounds .+ path_and_bounds[2]
    end
    return final_bounds
end


function handle_extra_edges_exact!(query::QueryGraph, data::DataGraph, partial_paths, current_query_nodes, visited_query_edges)
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
#            if query_edge_label == -1
#                weight *= length(data_edge_labels)
#            end
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
function get_exact_size(query::QueryGraph, data::DataGraph; use_partial_sums = true, verbose=false)
    node_order = topological_sort_by_dfs(bfs_tree(Graph(query.graph), vertices(query.graph)[1]))
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
    for node in vertices(data.graph)
        # if the id labels don't match, then don't initialize with this node
        query_data_label = get_data_label(query, new_node)
        if (query_data_label != -1)
            if (query_data_label != get_data_label(data, node))
                continue
            end
        end
        node_labels = data.vertex_labels[node]
        # if the node labels don't match, then don't initialize with this node
        if (parent_label == -1) || (in(parent_label, node_labels))
            push!(partial_paths, ([node], 1))
        end
    end
    while length(node_order) > 0
        handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges)
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
        query_child_id_label = query.vertex_id_labels[new_node]
        
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
                    if (query_child_id_label != -1)
                        if (query_child_id_label != data_child_id_label)
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
            else
                for data_new_node in inneighbors(data.graph, old_node)
                    new_weight = weight
                    # Only add a new partial path if the edge label and node label match our query.
                    data_edge_labels = data.edge_labels[(data_new_node,old_node)]
                    data_child_labels = data.vertex_labels[data_new_node]
                    data_child_id_label = get_data_label(data, data_new_node)
                    if (query_child_id_label != -1)
                        if (query_child_id_label != data_child_id_label)
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
    if verbose
        println("Current Query Nodes: ", current_query_nodes)
        println("Visited Query Edges: ", visited_query_edges)
        println("Number of Partial Paths: ", length(keys(partial_paths)))
    end
    handle_extra_edges_exact!(query, data, partial_paths, current_query_nodes, visited_query_edges)
    if verbose
        println("Current Query Nodes: ", current_query_nodes)
        println("Visited Query Edges: ", visited_query_edges)
        println("Number of Partial Paths: ", length(keys(partial_paths)))
    end

    
    exact_size = 0
    for path_and_weight in partial_paths 
        exact_size +=  path_and_weight[2]
    end 
    return exact_size
end
