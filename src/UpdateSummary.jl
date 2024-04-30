# This file contains a prototype implementation of updates to the lifted color summary of a Data Graph.

"""
Chooses a color from the lifted summary of the Data Graph.
Currently picks the color with the most nodes, but in future direction
different colors could be selected for node updates.
# Arguments
- summary::ColorSummary - the lifted summary describing the Data Graph.
"""
function choose_color(summary::ColorSummary)
    return get_largest_color(summary)
    # return get_update_only_color!(summary)
end

"""
Returns the color with the most nodes in the lifted Data Graph summary.
# Arguments
- summary::ColorSummary - the lifted summary describing the Data Graph.
"""
function get_largest_color(summary::ColorSummary)
    max_color_cardinality = 0
    current_color = collect(keys(summary.color_label_cardinality))[1] # initialize with the first color
    for color in keys(summary.color_label_cardinality)
        if summary.color_label_cardinality[color][-1] >= max_color_cardinality
            max_color_cardinality = summary.color_label_cardinality[color][-1]
            current_color = color
        end
    end
    return current_color
end

"""
"Adds" a node to the Data Graph by updating the summary statistics corresponding to the labels and colors of the added node.
# Arguments
- summary::ColorSummary{AvgDegStats} - the lifted summary describing the Data Graph.
- node_labels - the labels belonging to the node to add.
- node - the vertex ID of the node to add.
"""
function add_summary_node!(summary::ColorSummary{AvgDegStats}, node_labels, node)
    data_label = node - 1

    color = choose_color(summary)
    # add to the bloom filter
    push!(summary.color_filters[color], data_label)
    if Probably.loadfactor(summary.color_filters[color]) > .95
        push!(summary.color_full, color)
    end
    # for edge degrees, it decreases the average.
    # we want to update all the avg out/in degrees where this is the starting node (c1 == color && v1 == label),
    for edge_label in keys(summary.edge_deg)
        for node_label in node_labels
            if !haskey(summary.edge_deg[edge_label], node_label)
                summary.edge_deg[edge_label][node_label] = Dict()
            end
            if !haskey(summary.edge_deg[edge_label][node_label], color)
                summary.edge_deg[edge_label][node_label][color] = Dict()
            end
            for other_color in keys(summary.edge_deg[edge_label][node_label][color])
                current_ds = get(summary.edge_deg[edge_label][node_label][color], other_color, AvgDegStats(0, 0))
                current_cardinality = get(summary.color_label_cardinality[color], node_label, 0)
                avg_in = current_ds.avg_in * (current_cardinality / (current_cardinality + 1))
                avg_out = current_ds.avg_out * (current_cardinality / (current_cardinality + 1))
                new_ds = AvgDegStats(avg_in, avg_out)
                summary.edge_deg[edge_label][node_label][color][other_color] = new_ds
            end
        end
    end

    # add to the cardinality counts
    for node_label in node_labels
        summary.color_label_cardinality[color][node_label] = get(summary.color_label_cardinality[color], node_label, 0) + 1
    end
    if !(-1 in node_labels)
        summary.color_label_cardinality[color][-1] = get(summary.color_label_cardinality[color], -1, 0) + 1
    end
    summary.total_nodes += 1
end

# assumes that all nodes are currently in the color summary - the node is guaranteed to be in at least one bloom filter
"""
Finds and returns the color of the node in the lifted graph summary.
# Arguments
- summary::ColorSummary - the lifted color summary describing the Data Graph.
- node - the vertex ID of the node.
"""
function get_node_summary_color(summary::ColorSummary, node)
    possible_colors = collect(summary.color_full)
    for color in keys(summary.color_filters)
        filter = summary.color_filters[color]
        # in the data graph, the node's data label is just its id - 1
        if (node - 1) in filter
            push!(possible_colors, color)
        end
    end

    # Since Cuckoo filters are used, there is a chance that there will be false positive results.
    # In that case, we select the largest color.
    true_color = 0
    if length(possible_colors) == 0
        println("HERE")
        true_color = get_largest_color(summary)
    elseif length(possible_colors) == 1
        true_color = possible_colors[1]
    elseif length(possible_colors) > 1
        min_size = Inf
        for color in possible_colors
            color_size = summary.color_label_cardinality[color][-1]
            if color_size < min_size
                true_color = color
                min_size = color_size
            end
        end
    end

    return true_color
end

"""
"Adds" an edge to the Data Graph by updating the summary statistics corresponding to the edge labels and colors of the edge nodes.
# Arguments
- summary::ColorSummary{AvgDegStats} - the lifted Color Summary describing the Data Graph.
- start_node - the start node of the edge to add.
- end_node - the destination node of the edge to add.
- edge_labels::Vector - the labels belonging to the edge to add.
"""
function add_summary_edge!(summary::ColorSummary{AvgDegStats}, start_node, end_node, edge_labels::Vector)
    if !(-1 in edge_labels)
        push!(edge_labels, -1)
    end

    # first, find the distribution of labels for the current color pair
    start_color = get_node_summary_color(summary, start_node)
    end_color = get_node_summary_color(summary, end_node)
    for edge_label in edge_labels
        # will have to iterate over indegree vertex labels too, not just outdegree
        c1_count = summary.color_label_cardinality[start_color][-1]
        c2_count = summary.color_label_cardinality[end_color][-1]
        for vertex_label in keys(summary.color_label_cardinality[end_color])
            # this is the probability that the vertex label will be included
            probability_end_vertex_label = summary.color_label_cardinality[end_color][vertex_label] / summary.color_label_cardinality[end_color][-1]
            # the avg is based on the # of vertices in c1
            if !haskey(summary.edge_deg, edge_label)
                summary.edge_deg[edge_label] = Dict()
            end
            if !haskey(summary.edge_deg[edge_label], vertex_label)
                summary.edge_deg[edge_label][vertex_label] = Dict()
            end
            if !haskey(summary.edge_deg[edge_label][vertex_label], start_color)
                summary.edge_deg[edge_label][vertex_label][start_color] = Dict()
            end
            current_deg = get(summary.edge_deg[edge_label][vertex_label][start_color], end_color, AvgDegStats(0,0))
            original_avg_out = current_deg.avg_out
            new_avg_out = c1_count == 0 ? 0 :
            min(((original_avg_out * c1_count) + probability_end_vertex_label), c1_count * summary.color_label_cardinality[end_color][vertex_label]) / c1_count
            summary.edge_deg[edge_label][vertex_label][start_color][end_color] = AvgDegStats(current_deg.avg_in, new_avg_out)
            # note we don't have to update the color_label_cardinality since no new nodes were added...
        end
        for vertex_label in keys(summary.color_label_cardinality[start_color])
            probability_start_vertex_label = summary.color_label_cardinality[start_color][vertex_label] / summary.color_label_cardinality[start_color][-1]
            # now adjust the averages of the appropriate in/out degreestats...
            original_avg_in = 0
            if haskey(summary.edge_deg, edge_label) && haskey(summary.edge_deg[edge_label], vertex_label) &&
                    haskey(summary.edge_deg[edge_label][vertex_label], end_color) &&
                    haskey(summary.edge_deg[edge_label][vertex_label][end_color], start_color)
                original_avg_in = summary.edge_deg[edge_label][vertex_label][end_color][start_color].avg_in
            end
            if !haskey(summary.edge_deg, edge_label)
                summary.edge_deg[edge_label] = Dict()
            end
            if !haskey(summary.edge_deg[edge_label], vertex_label)
                summary.edge_deg[edge_label][vertex_label] = Dict()
            end
            if !haskey(summary.edge_deg[edge_label][vertex_label], end_color)
                summary.edge_deg[edge_label][vertex_label][end_color] = Dict()
            end
            current_deg = get(summary.edge_deg[edge_label][vertex_label][end_color], start_color, AvgDegStats(0, 0))
            original_avg_in = current_deg.avg_in
            new_avg_in =  c2_count == 0 ? 0 :
            min(((original_avg_in * c2_count) + probability_start_vertex_label), c2_count * summary.color_label_cardinality[start_color][vertex_label]) / c2_count
            summary.edge_deg[edge_label][vertex_label][end_color][start_color] = AvgDegStats(new_avg_in, current_deg.avg_out)
            # note we don't have to update the color_label_cardinality since no new nodes were added...
        end
    end
    summary.total_added_edges += 1
end
