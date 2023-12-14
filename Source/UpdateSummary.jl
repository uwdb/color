
# chooses a color for a new node to be added to
function choose_color(summary)
    # current implementation: find the biggest color
    # other future options:
    # - make a brand new color just for added nodes?
    return get_largest_color(summary)
end

function get_largest_color(summary)
    # current implementation: find the biggest color
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

function add_summary_node!(summary::ColorSummary{AvgDegStats}, node_labels, node)
    data_label = node 

    color = choose_color(summary)
    # add to the bloom filter
    push!(summary.color_filters[color], data_label)
    # for edge degrees, it decreases the average.
    # we want to update all the avg out degrees where this is the landing node (c2 == color && v2 == label),
    # and update all the avg in degrees where this is the starting node (c2 == color && v2 == label)
    for edge_label in keys(summary.edge_deg)
        for node_label in node_labels
            if !haskey(summary.edge_deg[edge_label], node_label)
                summary.edge_deg[edge_label][node_label] = Dict()
            end
            for other_color in keys(summary.edge_deg[edge_label][node_label])
                current_ds = get(summary.edge_deg[edge_label][node_label][other_color], color, AvgDegStats(0, 0))
                current_cardinality = get(summary.color_label_cardinality[color], node_label, 0)
                avg_in = current_ds.avg_in * (current_cardinality / (current_cardinality + 1))
                avg_out = current_ds.avg_out * (current_cardinality / (current_cardinality + 1))
                new_ds = AvgDegStats(avg_in, avg_out)
                summary.edge_deg[edge_label][node_label][other_color][color] = new_ds
            end
        end
    end

    # add to the cardinality counts
    for node_label in node_labels
        summary.color_label_cardinality[color][node_label] = get(summary.color_label_cardinality[color], node_label, 0) + 1
    end

    # for cycle stats, since the number of edges/cycles are the same,
    # cycle likelihood for an arbitrary edge doesn't change
end

# assume that you delete all attached edges before removing a summary node
# assume that the node to delete actually exists
function delete_summary_node!(summary::ColorSummary{AvgDegStats}, node_labels, node)
    color = get_node_summary_color(summary, node)

    # by definition in data graph
    data_label = node

    # remove from the cuckoo filter
    pop!(summary.color_filters[color], data_label)
    if (data_label in summary.color_filters[color])
        print("deletion from filter didn't work")
    end

    # adjust avg edge degrees
    for edge_label in keys(summary.edge_deg)
        for node_label in node_labels
            for other_color in keys(summary.edge_deg[edge_label][node_label])
                current_cardinality = get(summary.color_label_cardinality[color], node_label, 0)
                current_deg = get(summary.edge_deg[edge_label][node_label][other_color], color, AvgDegStats(0, 0))
                scale_factor = current_cardinality <= 1 ? 0 : (current_cardinality / (current_cardinality - 1))
                summary.edge_deg[edge_label][node_label][other_color][color] = AvgDegStats(current_deg.avg_in*scale_factor, current_deg.avg_out*scale_factor)
            end
        end
    end

    # subtract from the cardinality counts
    for node_label in node_labels
        summary.color_label_cardinality[color][node_label] -= 1
    end
end

# assumes that all nodes are currently in the color summary - the node is guaranteed to be in at least one bloom filter
function get_node_summary_color(summary, node)
    possible_colors = []
    for color in keys(summary.color_filters)
        filter = summary.color_filters[color]
        # in the data graph, the node's data label is just its id - 1
        if (node - 1) in filter
            push!(possible_colors, color)
        end
    end
    return length(possible_colors) == 0 ? rand(keys(summary.color_filters)) : rand(possible_colors)
end

function add_summary_edge!(summary, start_node, end_node, edge_labels)
    update_edge_degrees!(summary, start_node, end_node, edge_labels, remove=false)
end

function remove_summary_edge!(summary, start_node, end_node, edge_labels)
    update_edge_degrees!(summary, start_node, end_node, edge_labels, remove=true)
end

function update_edge_degrees!(summary::ColorSummary{AvgDegStats}, start_node, end_node, edge_labels::Vector; remove=false)
    # need to eventually make sure that the edge label is a set that includes -1 for the label...
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
            if remove
                probability_end_vertex_label *= -1
            end
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
            if remove
                probability_start_vertex_label *= -1
            end
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
