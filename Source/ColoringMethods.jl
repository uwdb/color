
function _quasi_stable_coloring(g::DataGraph, params::ColorSummaryParams, num_colors::Int)
    QSC = QuasiStableColors
    C = QSC.q_color(g.graph, n_colors=num_colors, weighting=params.weighting)
    color_hash::Dict{Int, Int} = QSC.node_map(C)
    return color_hash
end

function _hash_coloring(g::DataGraph, num_colors::Int)
    color_hash::Dict{Int, Int} = Dict()
    for i in 1:nv(g.graph)
        color_hash[i] = (hash(i) % num_colors) + 1
    end
    return color_hash
end

function _degree_coloring(g::DataGraph, num_colors::Int)
    color_hash = Dict()
    degrees = sort(degree(g.graph))
    bucket_right_edges = []
    for i in 1:num_colors
        degree_quantile = degrees[Int(floor(float(i)/num_colors *length(degrees)))]
        if i > 1 && bucket_right_edges[i-1] >= degree_quantile
            push!(bucket_right_edges, bucket_right_edges[i-1] + 1)
        else
            push!(bucket_right_edges, degree_quantile)
        end
    end
    for i in 1:nv(g.graph)
        node_degree = degree(g.graph, i)
        for j in 1:num_colors
            if node_degree <= bucket_right_edges[j]
                color_hash[i] = j
                break
            end
        end
    end
    return color_hash
end

function _directed_degree_coloring(g::DataGraph, num_colors::Int)
    color_hash = Dict()
    num_degree_buckets = Int(floor(float(num_colors)^.5))
    indegrees = sort(indegree(g.graph))
    in_bucket_right_edges = []
    for i in 1:num_degree_buckets
        degree_quantile = indegrees[Int(floor(float(i)/num_degree_buckets *length(indegrees)))]
        if i > 1 && in_bucket_right_edges[i-1] >= degree_quantile
            push!(in_bucket_right_edges, in_bucket_right_edges[i-1] + 1)
        else
            push!(in_bucket_right_edges, degree_quantile)
        end
    end

    outdegrees = sort(outdegree(g.graph))
    out_bucket_right_edges = []
    for i in 1:num_degree_buckets
        degree_quantile = outdegrees[Int(floor(float(i)/num_degree_buckets *length(outdegrees)))]
        if i > 1 && out_bucket_right_edges[i-1] >= degree_quantile
            push!(out_bucket_right_edges, out_bucket_right_edges[i-1] + 1)
        else
            push!(out_bucket_right_edges, degree_quantile)
        end
    end
    for i in 1:nv(g.graph)
        in_degree = indegree(g.graph, i)
        in_bucket = 0
        for j in 1:num_degree_buckets
            if in_degree <= in_bucket_right_edges[j]
                in_bucket = j
                break
            end
        end
        out_degree = outdegree(g.graph, i)
        out_bucket = 0
        for j in 1:num_degree_buckets
            if out_degree <= out_bucket_right_edges[j]
                out_bucket = j
                break
            end
        end
        color_hash[i] = (in_bucket - 1) * num_degree_buckets + out_bucket
    end
    return color_hash
end

# Group by labels, then combine groups of labels as needed
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


# another function: find the top n - 1 vertices with the most neighbors
# each of their neighbors are in that color, the rest are in a different color
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


# another function: color by label, sort labels by their avg number of edges, hash partition based on # edges
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

# This function takes in a coloring and attempts to refine each color into sub colors.
# It does this by recursively choosing a single label which most evenly breaks the color
# into two sub-colors with up to `label_refining_rounds` depth.
function _refine_by_vertex_labels(g::DataGraph, params::ColorSummaryParams,
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



# This function takes in a coloring and attempts to refine each color into sub colors.
# It does this by recursively choosing a single label which has the greatest stddev w.r.t.
# the edge count of vertices and splitting the nodes based on their edge count for that label
# with up to `label_refining_rounds` depth.
function _recursive_neighbor_split(g::DataGraph, group::Vector{NodeId}, depth::Int, max_depth::Int)
    if depth == max_depth
        return [group]
    end

    label_count_sequences = Dict{Int, Any}()
    for v in group
        label_counts = counter(Int)
        for n in all_neighbors(g.graph, v)
            for label in g.vertex_labels[n]
                inc!(label_counts, label)
            end
        end
        for label in keys(label_counts)
            if !haskey(label_count_sequences, label)
                label_count_sequences[label] = counter(Int)
            end
            inc!(label_count_sequences[label], label_counts[label])
        end
    end

    discriminating_label = -1
    discriminating_count = -1
    max_stddev = 0
    for label in keys(label_count_sequences)
        count_sequence = label_count_sequences[label]
        label_avg = sum([count_sequence[d]*d for d in keys(count_sequence)])/length(group)
        label_stddev = sum([(count_sequence[d]-label_avg)^2 for d in keys(count_sequence)])/length(group)
        label_stddev = sqrt(label_stddev)
        if label_stddev > max_stddev
            discriminating_label = label
            max_stddev = label_stddev
            discriminating_count = label_avg
        end
    end

    left_group = Vector{NodeId}()
    right_group = Vector{NodeId}()
    for v in group
        label_count = 0
        for n in all_neighbors(g.graph, v)
            if discriminating_label in g.vertex_labels[n]
                label_count += 1
            end
        end
        if label_count > discriminating_count
            push!(right_group, v)
        else
            push!(left_group, v)
        end
    end
    return [_recursive_neighbor_split(g, left_group, depth + 1, max_depth)..., _recursive_neighbor_split(g, right_group, depth + 1, max_depth)...]
end

function _refine_by_neighbor_labels(g::DataGraph, params::ColorSummaryParams,
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
        new_label_groups = _recursive_neighbor_split(g, color_to_vertices[c], 0, label_refining_rounds)
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

function color_graph(g::DataGraph, params::ColorSummaryParams, num_colors::Int)
    color_hash::Dict{NodeId, Color} = if params.partitioner == QuasiStable
         _quasi_stable_coloring(g, params, num_colors)
    elseif params.partitioner == Hash
         _hash_coloring(g, num_colors)
    elseif params.partitioner == Degree
         _degree_coloring(g, num_colors)
    elseif params.partitioner == DirectedDegree
         _directed_degree_coloring(g, num_colors)
    elseif params.partitioner == SimpleLabel
         _simple_label_coloring(g, num_colors)
    elseif params.partitioner == InOut
         _edge_ratio_color(g, num_colors)
    elseif params.partitioner == LabelInOut
         _label_edge_ratio_coloring(g, num_colors)
    elseif params.partitioner == NeighborEdges
         _label_edges_coloring(g, num_colors)
    elseif params.partitioner == MostNeighbors
         _label_most_neighbors_coloring(g, num_colors)
    end
    if params.label_refining_rounds > 0
        color_hash = _refine_by_neighbor_labels(g, params, color_hash, params.label_refining_rounds)
    end
    return color_hash
end
