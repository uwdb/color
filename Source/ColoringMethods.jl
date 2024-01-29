"""
Uses the Quasi-Stable coloring algorithm to (further) color a graph, partitioning nodes
based on their connections to nodes of other colors. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- params::ColorSummaryParams - parameters used for the Quasi-Stable coloring.
- color_hash::Dict{NodeId, Color} - pre-existing colors for the data graph, if any.
- new_colors::Int - the maximum number of new colors to add via Quasi-Stable coloring.
"""
function _quasi_stable_coloring(g::DataGraph, params::ColorSummaryParams, color_hash::Dict{NodeId, Color}, new_colors::Int)
    QSC = QuasiStableColors
    existing_colors = maximum(values(color_hash))
    color_vec = [NodeId[] for _ in 1:existing_colors]
    for (node, color) in color_hash
        push!(color_vec[color], node)
    end
    non_empty_vecs = [length(node_set) > 0 for node_set in color_vec]
    color_vec = color_vec[non_empty_vecs]
    C = QSC.q_color(g.graph, n_colors=existing_colors + new_colors, weighting=params.weighting, warm_start=color_vec)
    color_hash::Dict{Int, Int} = QSC.node_map(C)
    return color_hash
end

"""
Hash partitions data nodes by assigning colors at random. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- params::ColorSummaryParams - parameters used for the hash coloring.
- color_hash::Dict{NodeId, Color} - pre-existing colors for the data graph, if any.
- new_colors::Int - the maximum number of new colors to add via hash coloring.
"""
function _hash_coloring(g::DataGraph, params::ColorSummaryParams, color_hash::Dict{NodeId, Color}, new_colors::Int)
    existing_colors = maximum(values(color_hash))
    color_to_nodes = Dict{Color, Vector{NodeId}}(i => NodeId[] for i in 1:existing_colors)
    for (node, color) in color_hash
        push!(color_to_nodes[color], node)
    end
    split_color = 1
    next_color = existing_colors + 1
    colors_this_round = existing_colors
    while next_color < new_colors + existing_colors
        color_to_nodes[next_color] = []
        nodes_to_remove = copy(color_to_nodes[split_color])
        color_to_nodes[split_color] = []
        for node in nodes_to_remove
            if rand() < .5
                push!(color_to_nodes[next_color], node)
            else
                push!(color_to_nodes[split_color], node)
            end
        end

        if split_color == colors_this_round
            colors_this_round = next_color
            split_color = 1
        else
            split_color += 1
        end
        next_color += 1
    end
    color_hash = Dict()
    for (color, nodes) in color_to_nodes
        for node in nodes
            color_hash[node] = color
        end
    end
    return color_hash
end

"""
Partitions nodes based on their degrees. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- params::ColorSummaryParams - parameters used for the degree coloring.
- color_hash::Dict{NodeId, Color} - pre-existing colors for the data graph, if any.
- new_colors::Int - the maximum number of new colors to add via degree coloring.
"""
function _degree_coloring(g::DataGraph, params::ColorSummaryParams, color_hash::Dict{NodeId, Color}, new_colors::Int)
    existing_colors = maximum(values(color_hash))
    color_to_nodes = Dict{Color, Vector{NodeId}}(i => NodeId[] for i in 1:existing_colors)
    for (node, color) in color_hash
        push!(color_to_nodes[color], node)
    end

    function calc_max_dev_mean(c)
        degrees = Float64[]
        for node in color_to_nodes[c]
            push!(degrees, degree(g.graph, node))
        end
        deg_mean = mean(degrees)
        deg_total_dev = sum([(x-deg_mean)^2 for x in degrees])
        mean_deg = mean(degrees)
        return deg_total_dev, mean_deg
    end

    color_max_dev_mean = Dict{Color, Tuple{Float64, Float64}}()
    for c in 1:existing_colors
        color_max_dev_mean[c] = calc_max_dev_mean(c)
    end

    next_color = existing_colors + 1
    while next_color < new_colors + existing_colors
        split_color = 1
        split_mean_deg = -1
        max_deg_total_dev = -1
        for c in 1 : next_color - 1
            deg_total_dev, mean_deg = color_max_dev_mean[c]
            if deg_total_dev > max_deg_total_dev
                split_color = c
                max_deg_total_dev = deg_total_dev
                split_mean_deg = mean_deg
            end
        end
        nodes_to_remove = copy(color_to_nodes[split_color])
        color_to_nodes[next_color] = []
        color_to_nodes[split_color] = []
        for node in nodes_to_remove
            if degree(g.graph, node) <= split_mean_deg
                push!(color_to_nodes[next_color], node)
            else
                push!(color_to_nodes[split_color], node)
            end
        end
        if length(color_to_nodes[next_color]) > 0
            color_max_dev_mean[split_color] = calc_max_dev_mean(split_color)
            color_max_dev_mean[next_color] = calc_max_dev_mean(next_color)
            next_color += 1
        else
            break
        end
    end
    color_hash = Dict()
    for (color, nodes) in color_to_nodes
        for node in nodes
            color_hash[node] = color
        end
    end
    return color_hash
end

"""
Partitions nodes by grouping by label then combining groups until reaching the
limit for the number of new colors to add. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- num_colors::Int - the maximum number of new colors to add.
"""
function _simple_label_coloring(g::DataGraph, num_colors::Int)
    # returns a color hash mapping a node to its color
    color_hash = Dict()

    # first get a list of labels by iterating through all of the vertices
    overall_labels::Set{Int} = Set()
    for v in 1:nv(g.graph)
        current_labels = g.vertex_labels[v]
        for label in current_labels
            push!(overall_labels, label)
        end
    end
    label_groups::Vector{Vector{Int}} = [[x] for x in overall_labels]

    # condense the labels into groups until the number of "colors" is less than the requirement
    new_label_groups::Vector{Vector{Int}} = []
    num_groups = length(label_groups)
    while (num_groups > num_colors && num_groups >= 2)
        if length(label_groups) <= 1
            for label_group in new_label_groups
                push!(label_groups, label_group)
            end
            new_label_groups = []
        end
        first_label_group = pop!(label_groups)
        second_label_group = pop!(label_groups)
        # combine the two label groups we popped and push them to the new label groups
        for label in second_label_group
            push!(first_label_group, label)
        end
        push!(new_label_groups, first_label_group)
        num_groups = length(label_groups) + length(new_label_groups)
    end
    for label_group in new_label_groups
        push!(label_groups, label_group)
    end

    # for each label, map to an arbitrary color
    # label_colors[label] = color
    label_colors::Dict{Int, Int} = Dict()
    for x in eachindex(label_groups)
        for label in label_groups[x]
            label_colors[label] = x;
        end
    end

    # for each node in the graph, map it to its color based on an arbitrary label
    color_hash::Dict{Int, Int} = Dict()
    for v in 1:nv(g.graph)
        color_hash[v] = label_colors[g.vertex_labels[v][1]]
    end

    # return the color hash
    return color_hash
end

"""
Partitions nodes by assigning the top vertices with the most neighbors their own colors,
assigning their neighbors the corresponding colors, and putting any remaining nodes
in a separate color. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- num_colors::Int - the maximum number of new colors to add.
"""
function _label_most_neighbors_coloring(g::DataGraph, num_colors::Int)
    # map each node to how many vertices they are connected to
    vertex_mapping::Dict{Int, Int} = Dict()
    for v in 1:nv(g.graph)
        vertex_mapping[v] = length(neighbors(g.graph, v))
    end

    # sort into an array of node-neighbors pairs
    sorted_nodes = sort(collect(vertex_mapping), by=x->x[2], rev=true)

    # keep track of the top n - 1 vertices
    num_groups::Int = num_colors - 1

    color_hash = Dict()
    # for each node, go through the top n - 1 vertices to determine which group it belongs in (1...n -1)
    num_processed = 0
    for node_number_pair in vertex_mapping
        # if it doesn't neighbor any of the top vertices, its color will now be the "last" possible color
        current_color::Int = num_colors
        found_supernode = false
        for supernode in 1:num_groups
            if (!found_supernode)
                if (node_number_pair[1] == supernode) || (node_number_pair[1] in neighbors(g.graph, sorted_nodes[supernode][1]))
                    current_color = supernode
                end
            end
        end
        color_hash[node_number_pair[1]] = current_color
        num_processed+=1
    end
    return color_hash
end

"""
Partitions nodes by assigning arbitrary colors to specific node labels, then further partitioning
based on the number of connected edges. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- num_colors::Int - the maximum number of colors to add.
"""
function _label_edges_coloring(g::DataGraph, num_colors::Int)
    # returns a color hash mapping a node to its color
    color_hash = Dict()

    # map each label to the set of all the nodes with that label
    label_node_mapping::Dict{Int, Vector{Int}} = Dict()
    for v in 1:nv(g.graph)
        current_labels = g.vertex_labels[v]
        for label in current_labels
            if !(label in keys(label_node_mapping))
                label_node_mapping[label] = []
            end
            push!(label_node_mapping[label], v)
        end
    end

    # create a new mapping of label -> avg number of edges for each node within it
    label_edges_mapping::Dict{Int, Float64} = Dict()
    for label in keys(label_node_mapping)
        total_edges = 0
        for node in label_node_mapping[label]
            total_edges += length(neighbors(g.graph, node))
        end
        label_edges_mapping[label] = total_edges / length(label_node_mapping[label])
    end

    # now create a sorted list of labels based on their average number of edges
    sorted_labels = [x[1] for x in sort(collect(label_edges_mapping), by=x->x[2])]
    # now do very simple hash partition on the sorted list
    group_size::Int = (num_colors > length(sorted_labels)) ? 1 : ceil(Int, length(sorted_labels) / num_colors)

    # now for each node, assign a color based on the label group it belongs to
    for node in 1:nv(g.graph)
        labels = g.vertex_labels[node]
        # choose the label with the highest number of edges
        current_label = labels[1]
        for label in labels
            if label_edges_mapping[label] > label_edges_mapping[current_label]
                current_label = label
            end
        end
        color_hash[node] = (ceil(Int, indexin(current_label, sorted_labels)[1] / group_size))
    end
    return color_hash
end

"""
Partitions nodes by sorting nodes by label, sorting again by the average
indegree/outdegree ratio, then partitioning based on the sorted groups.
Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- num_colors::Int - the maximum number of colors to add.
"""
# another function: color by label, sort labels by avg in/out ratio, hash partition based on # colors
function _label_edge_ratio_coloring(g::DataGraph, num_colors::Int)
    # returns a color hash mapping a node to its color
    color_hash = Dict()

    # map each label to the set of all the nodes with that label (each node is randomly assigned to one of its labels)
    label_node_mapping::Dict{Int, Vector{Int}} = Dict()
    for v in 1:nv(g.graph)
        current_labels = g.vertex_labels[v]
        for label in current_labels
            if !(label in keys(label_node_mapping))
                label_node_mapping[label] = []
            end
            push!(label_node_mapping[label], v)
        end
    end

    # create a new mapping of label -> avg in/out ratio for each node within it
    label_edges_mapping::Dict{Int, Float64} = Dict()
    for label in keys(label_node_mapping)
        total_ratio::Float64 = 0
        for node in label_node_mapping[label]
            in_out_ratio::Float64 = length(inneighbors(g.graph, node)) / length(outneighbors(g.graph, node))
            total_ratio += in_out_ratio
        end
        label_edges_mapping[label] = total_ratio / length(label_node_mapping[label])
    end

    # now create a sorted list of labels based on their average in/out ratio
    sorted_labels = [x[1] for x in sort(collect(label_edges_mapping), by=x->x[2])]
    println("sorted labels")
    # now do very simple hash partition on the sorted list
    group_size::Int = (num_colors > length(sorted_labels)) ? 1 : ceil(Int, length(sorted_labels) / num_colors)

    # now for each node, assign a color based on the label group it belongs to
    for node in 1:nv(g.graph)
        labels = g.vertex_labels[node]
        # choose the label with the highest number of edges
        current_label = labels[1]
        for label in labels
            if label_edges_mapping[label] > label_edges_mapping[current_label]
                current_label = label
            end
        end
        color_hash[node] = (ceil(Int, indexin(current_label, sorted_labels)[1] / group_size))
    end

    # now return the color mapping.
    return color_hash
end

"""
Partitions nodes by sorting based on indegree/outdegree ratio then assigning
colors to groups with similar ratios. Returns a mapping of nodes to their colors.
# Arguments
- g::Datagraph - the data graph to color.
- num_colors::Int - the maximum number of colors to add.
"""
# edge ratio without grouping by labels beforehand
function _edge_ratio_color(g::DataGraph, num_colors::Int)
    # returns a color hash mapping a node to its color
    color_hash = Dict()

    # map each node to its average in/out ratio
    node_ratio_mapping::Dict{Int, Float64} = Dict()
    for v in 1:nv(g.graph)
        in_out_ratio::Float64 = length(inneighbors(g.graph, v)) / length(outneighbors(g.graph, v))
        node_ratio_mapping[v] = in_out_ratio
    end

    # sort the nodes by their in/out ratio
    sorted_nodes = [x[1] for x in sort(collect(node_ratio_mapping), by=x->x[2])]

    # create a color hash based on the sorted nodes and the number of desired colors
    group_size::Int = (num_colors > length(sorted_nodes)) ? 1 : ceil(Int, length(sorted_nodes) / num_colors)
    for v in 1:nv(g.graph)
        color_hash[v] = (ceil(Int, indexin(v, sorted_nodes)[1] / group_size))
    end
    # now return the color mapping.
    return color_hash
end

"""
Recursively partitions nodes based on which label is most evenly split between those
in the group that do (not) have the label. Returns a list of node groupings.
# Arguments
- g::Datagraph - the data graph to color
- group::Vector{NodeId} - the current group of nodes to num_process
- depth::Int - current depth of recursive processing
- max_depth::Int - maximum depth of recursive processing.
"""
function _recursive_label_split(g::DataGraph, group::Vector{NodeId}, depth::Int, max_depth::Int)
    if depth == max_depth
        return [group]
    end

    label_counts = counter(Int)
    for v in group
        for label in g.vertex_labels[v]
            inc!(label_counts, label)
        end
    end

    most_even_label = -1
    max_evenness = 1.0
    group_size = float(length(group))
    for l in keys(label_counts)
        if abs(label_counts[l]/group_size - .5) < max_evenness
            max_evenness = abs(label_counts[l]/group_size - .5)
            most_even_label = l
        end
    end

    left_group = Vector{NodeId}()
    right_group = Vector{NodeId}()
    for v in group
        if most_even_label in g.vertex_labels[v]
            push!(right_group, v)
        else
            push!(left_group, v)
        end
    end
    return [_recursive_label_split(g, left_group, depth + 1, max_depth)..., _recursive_label_split(g, right_group, depth + 1, max_depth)...]
end

"""
Takes in a coloring and attempts to refine each color into sub colors.
Does this by recursively choosing a single label which most evenly breaks the color
into two sub-colors.
# Arguments
- g::DataGraph - the data graph to color.
- color_hash::Dict{NodeId, Color} - the existing coloring to refine.
- label_refining_rounds::Int - the number of times to recursively refine the color.
"""
function _refine_by_vertex_labels(g::DataGraph,
                                    color_hash::Dict{NodeId, Color}, label_refining_rounds::Int)
    refined_color_hash::Dict{NodeId, Color} = Dict()
    color_to_vertices::Dict{Color, Vector{NodeId}} = Dict()
    for v in keys(color_hash)
        color = color_hash[v]
        if haskey(color_to_vertices, color)
            push!(color_to_vertices[color], v)
        else
            color_to_vertices[color] = [v]
        end
    end

    color_counter = 1
    for c in keys(color_to_vertices)
        new_label_groups = _recursive_label_split(g, color_to_vertices[c], 0, label_refining_rounds)
        for group in new_label_groups
            if length(group) == 0
                continue
            end
            for v in group
                refined_color_hash[v] = color_counter
            end
            color_counter += 1
        end
    end
    return refined_color_hash
end

"""
Takes in a coloring and attempts to refine it to add `new_colors` number of colors. 
Does this by choosing a single label which has the greatest stddev w.r.t.
the edge count of vertices and splitting the nodes based on their edge count for that label.
# Arguments
- g::DataGraph - the data graph to color.
- color_hash::Dict{NodeId, Color} - the existing coloring to refine.
- new_colors::Int - the number of new colors to add via splitting by edge count.
"""
function _neighbor_labels_coloring(g::DataGraph, params::ColorSummaryParams,
                                        color_hash::Dict{NodeId, Color}, new_colors::Int)
    existing_colors = maximum(values(color_hash))
    color_to_nodes = Dict{Color, Vector{NodeId}}(i => NodeId[] for i in 1:existing_colors)
    for (node, color) in color_hash
        push!(color_to_nodes[color], node)
    end
    next_color = existing_colors + 1

    node_label_counts = Dict()
    for node in keys(color_hash)
        node_label_counts[node] = counter(Int)
        for neighbor in all_neighbors(g.graph, node)
            for label in g.vertex_labels[neighbor]
                inc!(node_label_counts[node], label)
            end
        end
    end


    function calc_max_dev_label_mean(c)
        split_label = -1
        max_deg_total_dev = -1
        split_mean_deg = -1
        label_degrees = Dict()
        for node in color_to_nodes[c]
            label_counts = node_label_counts[node]
            for (label, count) in label_counts
                if !haskey(label_degrees, label)
                    label_degrees[label] = Int[]
                end
                push!(label_degrees[label], count)
            end
        end
        for (label, degrees) in label_degrees
            degrees = union(degrees, zeros(length(color_to_nodes[c])-length(degrees)))
            deg_mean = mean(degrees)
            deg_total_dev = sum([(x-deg_mean)^2 for x in degrees])
            if deg_total_dev > max_deg_total_dev
                split_label = label
                max_deg_total_dev = deg_total_dev
                split_mean_deg = mean(degrees)
            end
        end
        return max_deg_total_dev, split_label, split_mean_deg
    end

    color_max_dev_label_mean = Dict{Color, Tuple{Float64, Int, Float64}}()
    for c in 1:existing_colors
        color_max_dev_label_mean[c] = calc_max_dev_label_mean(c)
    end

    while next_color < new_colors + existing_colors
        split_color = 1
        split_label = 1
        split_mean_deg = -1
        max_deg_total_dev = -1
        for c in 1 : next_color - 1
            if haskey(color_max_dev_label_mean, c)
                deg_total_dev, label, mean_deg = color_max_dev_label_mean[c]
                if deg_total_dev > max_deg_total_dev
                    split_color = c
                    split_label = label
                    max_deg_total_dev = deg_total_dev
                    split_mean_deg = mean_deg
                end
                continue
            end
        end
        nodes_to_remove = copy(color_to_nodes[split_color])
        color_to_nodes[next_color] = []
        color_to_nodes[split_color] = []
        for node in nodes_to_remove
            node_degree = node_label_counts[node][split_label]
            if node_degree <= split_mean_deg
                push!(color_to_nodes[next_color], node)
            else
                push!(color_to_nodes[split_color], node)
            end
        end
        color_max_dev_label_mean[split_color] = calc_max_dev_label_mean(split_color)
        color_max_dev_label_mean[next_color] = calc_max_dev_label_mean(next_color)
        next_color += 1
    end

    color_hash = Dict()
    for (color, nodes) in color_to_nodes
        for node in nodes
            color_hash[node] = color
        end
    end
    return color_hash
end

"""
Takes in a coloring and attempts to refine it to add `new_colors` number of colors. 
Does this by choosing finding the label with the most even distribution of nodes that do (not)
have it then splitting based on this label.
# Arguments
- g::DataGraph - the data graph to color.
- color_hash::Dict{NodeId, Color} - the existing coloring to refine.
- new_colors::Int - the number of new colors to add via splitting by edge count.
"""
function _vertex_labels_coloring(g::DataGraph,
                                        color_hash::Dict{NodeId, Color}, new_colors::Int)
    existing_colors = maximum(values(color_hash))
    color_to_nodes = Dict{Color, Vector{NodeId}}(i => NodeId[] for i in 1:existing_colors)
    for (node, color) in color_hash
        push!(color_to_nodes[color], node)
    end
    next_color = existing_colors + 1

    function calc_most_avg_label_split(c)
        split_label = -1
        split_equity = 1.0
        num_nodes = length(color_to_nodes[c])
        label_counts = counter(Int)
        for node in color_to_nodes[c]
            for label in g.vertex_labels[node]
                inc!(label_counts, label)
            end
        end
        for (label, count) in label_counts
            label_prop = count / num_nodes
            if abs(label_prop - .5) < split_equity
                split_label = label
                split_equity = abs(label_prop - .5)
            end
        end
        return split_equity, split_label
    end

    color_most_avg_label_equity = Dict{Color, Tuple{Float64, Int}}()
    for c in 1:existing_colors
        color_most_avg_label_equity[c] = calc_most_avg_label_split(c)
    end

    while next_color < new_colors + existing_colors
        split_color = 1
        split_label = 1
        weighted_split_equity = Inf
        for c in 1 : next_color - 1
            if haskey(color_most_avg_label_equity, c)
                equity, label = color_most_avg_label_equity[c]
                weighted_equity = (.5 - equity) * length(color_to_nodes[c])
                if weighted_equity < weighted_split_equity
                    split_color = c
                    split_label = label
                    weighted_split_equity = weighted_equity
                end
            end
        end
        nodes_to_remove = copy(color_to_nodes[split_color])
        color_to_nodes[next_color] = []
        color_to_nodes[split_color] = []
        for node in nodes_to_remove
            has_label = split_label in g.vertex_labels[node]
            if has_label
                push!(color_to_nodes[next_color], node)
            else
                push!(color_to_nodes[split_color], node)
            end
        end
        color_most_avg_label_equity[split_color] = calc_most_avg_label_split(split_color)
        color_most_avg_label_equity[next_color] = calc_most_avg_label_split(next_color)
        next_color += 1
    end

    color_hash = Dict()
    for (color, nodes) in color_to_nodes
        for node in nodes
            color_hash[node] = color
        end
    end
    return color_hash
end

"""
Colors a data graph then continues refining the colors as described by the given parameters. Returns a mapping of nodes to their assigned color.
# Arguments
- g::DataGraph - the data graph to color.
- params::ColorSummaryParams - parameters describing how to color the graph.
# Example
    color_graph(g, new ColorSummaryParams(partioning_scheme=[(QuasiStable, 8), (Hash, 8)]))
The above method call will first do a 8-color Quasi-Stable coloring of a DataGraph g, then refine with 8 more colors via a hash partitioning.

"""
function color_graph(g::DataGraph, params::ColorSummaryParams)
    if nv(g.graph) == 0
        return Dict()
    end
    color_hash::Dict{NodeId, Color} = Dict(i => 1 for i in 1:nv(g.graph))
    for (partitioner, num_colors) in params.partitioning_scheme
        color_hash = if partitioner == QuasiStable
                _quasi_stable_coloring(g, params, color_hash, num_colors)
        elseif partitioner == Hash
                _hash_coloring(g, params, color_hash, num_colors)
        elseif partitioner == Degree
                _degree_coloring(g, params, color_hash, num_colors)
        elseif partitioner == NeighborNodeLabels
                _neighbor_labels_coloring(g, params, color_hash, num_colors)
        elseif partitioner == NodeLabels
                _vertex_labels_coloring(g, color_hash, num_colors)
#        elseif partitioner == SimpleLabel
#                _simple_label_coloring(g, params, color_hash, num_colors)
#        elseif partitioner == InOut
#                _edge_ratio_color(g, params, color_hash, num_colors)
#        elseif partitioner == LabelInOut
#                _label_edge_ratio_coloring(g, params, color_hash, num_colors)
#        elseif partitioner == NeighborEdges
#                _label_edges_coloring(g, params, color_hash, num_colors)
#        elseif partitioner == MostNeighbors
#                _label_most_neighbors_coloring(g, params, color_hash, num_colors)
        else
            throw(ErrorException(string(partitioner) * " Not Implemented"))
        end
    end
    return color_hash
end
