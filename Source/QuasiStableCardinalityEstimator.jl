# This file contains a prototype implementation of Quasi-Stable Cardinality Estimation.

# The following two functions sum over all paths which have the same color assigned to a particular node in the query graph.
# Equivalently, they perform a groupby on all other nodes of the query graph. The goal of this is to prevent
# an exponential growth in the number of paths through the lifted color graph. However, we can only remove query nodes whose
# edges have already been processed.
@inline function sum_over_node(partial_paths::Matrix{Color}, partial_weights::Matrix{Float64}, current_query_nodes, node_to_remove)
    nodeIdx = 1
    for node in current_query_nodes
        if node == node_to_remove
            break
        end
        nodeIdx += 1
    end
    new_partial_paths::Dict{Vector{Color}, Vector{Float64}} = Dict()
    for i in 1:size(partial_paths)[2]
        new_path = copy(partial_paths[:, i])
        deleteat!(new_path, nodeIdx)
        if !haskey(new_partial_paths, new_path)
            new_partial_paths[new_path] = partial_weights[:, i]
        else
            new_partial_paths[new_path] .+= partial_weights[:, i]
        end
    end
    deleteat!(current_query_nodes, nodeIdx)
    partial_paths = zeros(Color, length(current_query_nodes), length(keys(new_partial_paths)))
    partial_weights = zeros(Float64, 3, length(keys(new_partial_paths)))

    path_idx = 1
    for path in keys(new_partial_paths)
        for i in 1:length(path)
            partial_paths[i, path_idx] = path[i]
        end
        weights = new_partial_paths[path]
        partial_weights[1, path_idx] = weights[1]
        partial_weights[2, path_idx] = weights[2]
        partial_weights[3, path_idx] = weights[3]
        path_idx += 1
    end
    return partial_paths, partial_weights
end

@enum SAMPLING_STRATEGY uniform weighted redistributive online loop_vec

@inline function sample_paths(partial_paths::Matrix{Color}, partial_weights::Matrix{Float64}, num_samples::Int, sampling_strategy::SAMPLING_STRATEGY)
    # if we want to sample more paths than there are existing, then just return the original partial paths
    num_nonzero_entries = 0
    for i in 1:size(partial_weights)[2]
        if partial_weights[2, i] > 0
            num_nonzero_entries += 1
        end
    end
    if (num_samples > num_nonzero_entries) || sampling_strategy == online
        return partial_paths, partial_weights
    end


    # sum up all of the bounds
    overall_bounds_sum::Float64 = sum(partial_weights[2, :])
    # choose a sample of the paths
    sample_weights = partial_weights[2, :]
    sample_weights = AnalyticWeights(sample_weights ./ overall_bounds_sum)
    if sampling_strategy == uniform
        sample_weights = AnalyticWeights([partial_weights[2, i] > 0 ? 1.0 : 0.0 for i in 1:size(partial_paths)[2]] ./ num_nonzero_entries)
    end
    sample_indices::Vector{Int} = sample(1:size(partial_weights)[2], sample_weights,  num_samples; replace=false)

    # sum up the sampled bounds
    sampled_bounds_sum::Float64 = 0
    for idx in sample_indices
        sampled_bounds_sum += partial_weights[2, idx]
    end

    # get the difference between the overall and sampled bound sum_over_finished_query_nodes
    bound_diff::Float64 = overall_bounds_sum - sampled_bounds_sum

    new_partial_paths = zeros(Color, size(partial_paths)[1], length(sample_indices))
    new_partial_weights = zeros(Float64, 3, length(sample_indices))

    # for each sampled path...
    new_path_idx = 1
    for idx in sample_indices
        new_partial_paths[:, new_path_idx] .= partial_paths[:, idx]
        if sampling_strategy == redistributive || sampling_strategy == loop_vec
            # figure out what fraction of the sampled bounds is in the current Bounds
            # higher bounds will have more weight redistributed to them
            bound_fractions = partial_weights[2, idx] / sampled_bounds_sum
            # use that fraction of the difference (i.e. the removed path weights) and add it to the partial path
            redistributed_weights = bound_fractions * bound_diff
            new_partial_weights[:, new_path_idx] .= partial_weights[:, idx] .+ redistributed_weights
        elseif sampling_strategy == weighted
            inverse_sampling_probability = overall_bounds_sum /  partial_weights[2, idx] / length(sample_indices)
            new_partial_weights[:, new_path_idx] .= partial_weights[:, idx] .* inverse_sampling_probability
        elseif sampling_strategy == uniform
            inverse_sampling_probability = size(partial_paths)[2] / length(sample_indices)
            new_partial_weights[:, new_path_idx] .= partial_weights[:, idx] .* inverse_sampling_probability
        end
        new_path_idx += 1
    end
    return new_partial_paths, new_partial_weights
end

function get_simple_paths_dfs!(visited::Set{Int}, cur::Int, finish::Int, max_length::Int,
                                graph::SimpleGraph, current_path::Vector{Int},
                                simple_paths::Vector{Vector{Int}})
    length(current_path) > max_length && return
    cur in visited && return
    push!(visited, cur)
    push!(current_path, cur)
    if cur == finish
        push!(simple_paths, deepcopy(current_path))
        delete!(visited, cur)
        pop!(current_path)
        return
    end

    for next in all_neighbors(graph, cur)
        get_simple_paths_dfs!(visited, next, finish, max_length, graph, current_path,
                                simple_paths)
    end
    if length(current_path) > 0
        pop!(current_path)
    end
    delete!(visited, cur)
end

# gets all directed, simple paths from the start to finish node
function get_all_simple_path_bools(start::Int, finish::Int, max_length::Int,
                                    query_graph::DiGraph, visited_edges::Vector{Tuple{Int,Int}})
    # convert the graph to be undirected and only include the edges that have already been processed
    graph_copy = Graph(nv(query_graph))
    for edge in visited_edges
        add_edge!(graph_copy, edge[1], edge[2])
    end
    rem_edge!(graph_copy, start, finish)

    visited = Set{Int}()
    current_path = Vector{Int}()
    simple_paths = Vector{Vector{Int}}()
    get_simple_paths_dfs!(visited, start, finish, max_length, graph_copy,
                                current_path, simple_paths)
    path_bools = Vector{BoolPath}()
    for path in simple_paths
        bools::Vector{Bool} = [false for _ in 1:length(path)-1]
        for i in 1 : length(path)-1
            src_node = path[i]
            dst_node = path[i+1]
            if dst_node in outneighbors(query_graph, src_node)
                bools[i] = true # out edge
            else
                bools[i] = false # in edge
            end
        end
        push!(path_bools, bools)
    end
    return path_bools
end

# gets the directed path from the start to finish node
function get_matching_graph(start::Int, finish::Int, query::QueryGraph)
    # convert the graph to be undirected
    graph_copy = Graph(copy(query.graph))
    rem_edge!(graph_copy, start, finish)
    # get a path from the start to finish node
    edges = a_star(graph_copy, start, finish)
    # currently assumes the edges go in order of the path,
    # will have to debug through and see if it works
    new_graph = DiGraph(length(edges) + 1)
    current_start = 0
    for edge in edges
        current_start += 1
        if src(edge) in outneighbors(query.graph, dst(edge))
            # this is a backwards edge
            add_edge!(new_graph, current_start + 1, current_start)
        else
            # this is a forwards edge
            add_edge!(new_graph, current_start, current_start + 1)
        end
    end
    return new_graph
end

@inline function handle_extra_edges!(query::QueryGraph, summary::ColorSummary, partial_paths::Array{Color}, partial_weights::Array{Float64},
                                current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int,Int}}, usingStoredStats::Bool,
                                only_shortest_path_cycle::Bool)
    # To account for cyclic queries, we check whether there are any remaining edges that have not
    # been processed. If so, we set the lower bound to 0, reduce the average estimate accordingly, and leave
    # the upper bound unchanged.
    remaining_edges::Vector{Tuple{Int, Int}} = []
    for edge in edges(query.graph)
        if ! ((src(edge), dst(edge)) in visited_query_edges) &&
                 (src(edge) in current_query_nodes && dst(edge) in current_query_nodes)
            push!(remaining_edges, (src(edge), dst(edge)))
        end
    end

    # scale down the average if there are remaining non-tree-edges
    for edge in remaining_edges
        push!(visited_query_edges, edge)
        parent_node_idx::Int = only(indexin(edge[1], current_query_nodes))
        new_node_idx::Int = only(indexin(edge[2], current_query_nodes))
        child_label::Int = only(query.vertex_labels[edge[2]])
        edge_label::Int = only(query.edge_labels[(edge[1],edge[2])])
        all_path_bools = Vector{BoolPath}()
        if only_shortest_path_cycle
            all_path_bools = [convert_path_graph_to_bools(get_matching_graph(edge[2], edge[1], query))]
        else
            all_path_bools = get_all_simple_path_bools(edge[2], edge[1], summary.max_cycle_size, query.graph, visited_query_edges)
        end

        default_colors::StartEndColorPair = (-1, -1)
        edge_deg::Dict{Int, Dict{Int, DegreeStats}} = Dict()
        if haskey(summary.edge_deg, edge_label) && haskey(summary.edge_deg[edge_label], child_label)
            edge_deg = summary.edge_deg[edge_label][child_label]
        end
        for i  in 1:size(partial_paths)[2]
            parent_color::Color = partial_paths[parent_node_idx, i]
            child_color::Color =  partial_paths[new_node_idx, i]
            current_colors::StartEndColorPair = (child_color, parent_color)
            # We don't have to check data label because these nodes are already in the
            # partial path, so we have already ensured that the colors are appropriate
            probability_no_edge = 1.0
            if (haskey(edge_deg, parent_color) && haskey(edge_deg[parent_color], child_color))
                if usingStoredStats && length(all_path_bools) > 0
                    for path_bools in all_path_bools
                        path_length = length(path_bools)
                        default_cycle_description = CyclePathAndColors(path_bools, default_colors)
                        current_cycle_description = CyclePathAndColors(path_bools, current_colors)
                        if haskey(summary.cycle_probabilities, current_cycle_description)
                            probability_no_edge *= 1.0 - summary.cycle_probabilities[current_cycle_description]
                        elseif haskey(summary.cycle_probabilities, default_cycle_description)
                            probability_no_edge *= 1.0 - summary.cycle_probabilities[default_cycle_description]
                        elseif haskey(summary.cycle_length_probabilities, path_length)
                            probability_no_edge *= 1.0 - summary.cycle_length_probabilities[path_length]
                        else
                            probability_no_edge *= 1.0 - get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary)
                        end
                    end
                else
                    probability_no_edge *= 1.0 - get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary)
                end
            end
            partial_weights[1, i] = 0
            partial_weights[2, i] *=  (1.0 - probability_no_edge)
        end
    end
end

function sum_over_finished_query_nodes(query::QueryGraph, partial_paths::Matrix{Color}, partial_weights::Matrix{Float64},
                                            current_query_nodes::Vector{Int}, visited_query_edges::Vector{Tuple{Int, Int}})
    new_partial_paths, new_partial_weights = partial_paths, partial_weights
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
            new_partial_paths, new_partial_weights = sum_over_node(new_partial_paths, new_partial_weights, current_query_nodes, node)
        end
    end
    return new_partial_paths, new_partial_weights
end


function get_cardinality_bounds(query::QueryGraph, summary::ColorSummary; max_partial_paths::Union{Nothing, Int} = nothing,
                                use_partial_sums::Bool = true, verbose::Bool = false, usingStoredStats::Bool = false,
                                include_cycles::Bool = true, sampling_strategy::SAMPLING_STRATEGY=weighted,
                                only_shortest_path_cycle::Bool=false)
    node_order::Vector{Int} = get_min_width_node_order(query.graph) #spanning tree to cut out cycles
    if verbose
        println("Node Order:", node_order)
    end
    # Because the label is implied by the color -> query_graph_vertex mapping stored in current_query_nodes,
    # we don't have to keep the label in the partial paths object.
    num_colors = length(summary.color_label_cardinality)
    partial_paths = zeros(Color, 1, num_colors) # each tuple contains a pairing of color paths -> bounds
    partial_weights = zeros(Float64, 3, num_colors)
    visited_query_edges::Vector{Tuple{Int,Int}} = []
    current_query_nodes::Vector{Int} = []

    old_node = popfirst!(node_order)
    parent_label = only(query.vertex_labels[old_node])
    parent_data_labels = get_data_label(query, old_node)
    push!(current_query_nodes, old_node)
    # Initialize partial_paths with all possible starting color/vertex possibilities.
    for color in 1:num_colors
        partial_paths[1, color] = color
        !haskey(summary.color_filters, color) && continue
        # Only use the parent label.
        # if the parent has a specified data label, only use colors that the filters approve
        data_label_is_in_color = false
        for data_label in parent_data_labels
            if data_label == -1
                data_label_is_in_color = true
                continue
            end
            if data_label in summary.color_filters[color]
                data_label_is_in_color = true
            end
        end
        partial_paths[1, color] = color
        if data_label_is_in_color && haskey(summary.color_label_cardinality[color], parent_label)
            partial_weights[1, color] = summary.color_label_cardinality[color][parent_label]
            partial_weights[2, color] = summary.color_label_cardinality[color][parent_label]
            partial_weights[3, color] = summary.color_label_cardinality[color][parent_label]
        end
    end
    rng = Random.default_rng()
    new_node = old_node
    while length(node_order) > 0
        if verbose
            println("Current Query Nodes: ", current_query_nodes)
            println("Visited Query Edges: ", visited_query_edges)
            println("Number of Partial Paths: ", length(keys(partial_paths)))
        end
        if use_partial_sums
            partial_paths, partial_weights = sum_over_finished_query_nodes(query, partial_paths, partial_weights, current_query_nodes, visited_query_edges)
        end
        if verbose
            println("Number of Partial Paths After Sum: ", length(keys(partial_paths)))
            println("Current Query Nodes After Sum: ", current_query_nodes)
            println("Visited Query Edges After Sum: ", visited_query_edges)
        end

        # Get the next child from the node order.
        new_node = popfirst!(node_order)
        parent_idx = 0
        out_edge = false
        for neighbor in all_neighbors(query.graph, new_node)
            if neighbor in current_query_nodes
                old_node = neighbor
                parent_idx::Int = only(indexin(neighbor, current_query_nodes))
                if old_node in inneighbors(query.graph, new_node)
                    out_edge = true
                end
                break
            end
        end
        # Push the current edge and nodes to the visited lists.
        if out_edge
            push!(visited_query_edges, (old_node, new_node))
        else
            push!(visited_query_edges, (new_node, old_node))
        end
        push!(current_query_nodes, new_node)

        # Get the appropriate labels, the query only uses one label per vertex/node.
        edge_label = 0
        if out_edge
            edge_label = only(query.edge_labels[(old_node,new_node)])
        else
            edge_label =  only(query.edge_labels[(new_node,old_node)])
        end

        new_label = only(query.vertex_labels[new_node])
        new_data_labels::Vector{Int} = get_data_label(query, new_node)
        edge_deg::Dict{Color, Dict{Color, DegreeStats}} = Dict()
        if haskey(summary.edge_deg, edge_label) &&
                        haskey(summary.edge_deg[edge_label], new_label)
            edge_deg = summary.edge_deg[edge_label][new_label]
        end
        num_old_paths = size(partial_paths)[2]
        new_path_idx = 1
        total_weight = 0
        if sampling_strategy == online
            num_new_paths = max_partial_paths
            new_partial_paths = zeros(Color, length(current_query_nodes), num_new_paths)
            new_partial_weights = zeros(Float64, 3, num_new_paths)
            H = PriorityQueue{Tuple{Int, Color, Tuple{Float64,Float64,Float64}}, Float64}()
            # Update the partial paths using the parent-child combo that comes next from the query.
            X = 0.0
            for i in 1:num_old_paths
                old_color = partial_paths[parent_idx, i]
                # Account for colors with no outgoing children.
                if haskey(edge_deg, old_color)
                    for new_color in keys(edge_deg[old_color])
                        # revamp the logic to use a set of labels rather than just one
                        # check if the data label(s) are in the color
                        data_label_in_color = false
                        for data_label in new_data_labels
                            if data_label == -1 ||  data_label in summary.color_filters[new_color]
                                data_label_in_color = true
                                break
                            end
                        end
                        !data_label_in_color && continue

                        degree_stats::DegreeStats = edge_deg[old_color][new_color]
                        new_min::Float64 = partial_weights[1, i] * (out_edge ? degree_stats.min_out : degree_stats.min_in)
                        new_avg::Float64 = partial_weights[2, i] * (out_edge ? degree_stats.avg_out : degree_stats.avg_in)
                        new_max::Float64 = partial_weights[3, i] * (out_edge ? degree_stats.max_out : degree_stats.max_in)
                        if !(length(new_data_labels) == 1 && new_data_labels[1] == -1)
                            # we have already confirmed that the data label is in the color, but if the data label isn't -1
                            # then we need to scale down the result since we only want to consider one of the many nodes in the new color
                            # we also need to set the minimum to 0 but keep the maximum the same
                            new_min = 0
                            new_avg /= summary.color_label_cardinality[new_color][new_label]
                        end
                        total_weight += new_avg
                        if length(H) < num_new_paths
                            heap_val = rand(rng)^(1/new_avg)
                            enqueue!(H, (i, new_color, (new_min, new_avg, new_max))=> heap_val)
                            if length(H) == num_new_paths
                                X = log(rand(rng)) / log(first(H)[2])
                            end
                        else
                            X = X-new_avg
                            if X <= 0
                                t = first(H)[2] ^ new_avg
                                r = (rand(rng)*(1-t) + t) ^ (1/new_avg)
                                dequeue!(H)
                                enqueue!(H, (i, new_color, (new_min, new_avg, new_max)) => r)
                                X = log(rand(rng)) / log(first(H)[2])
                            end
                        end
                    end
                end
            end
            for path_color_weights in keys(H)
                path_idx = path_color_weights[1]
                new_color = path_color_weights[2]
                weights = path_color_weights[3]
                for j in 1:length(current_query_nodes)-1
                    new_partial_paths[j, new_path_idx] = partial_paths[j, path_idx]
                end
                new_partial_paths[length(current_query_nodes), new_path_idx] = new_color
                new_partial_weights[1, new_path_idx] = weights[1]
                new_partial_weights[2, new_path_idx] = weights[2]
                new_partial_weights[3, new_path_idx] = weights[3]
                new_path_idx += 1
            end
            sample_sum = sum(new_partial_weights[2, :])
            if sample_sum > 0
                new_partial_weights = new_partial_weights .* total_weight/sample_sum
            end
            partial_paths = new_partial_paths[:, 1:new_path_idx-1]
            partial_weights = new_partial_weights[:, 1:new_path_idx-1]
        elseif sampling_strategy == loop_vec
            new_partial_paths = zeros(Color, length(current_query_nodes),  num_colors, num_old_paths)
            new_partial_weights = zeros(Float64, 3,  num_colors, num_old_paths)
            new_degs = zeros(Float32, num_colors, num_colors, 3)
            for i in keys(edge_deg)
                for (j, stat) in edge_deg[i]
                    if out_edge
                        new_degs[i,j,1] = stat.min_out
                        new_degs[i,j,2] = stat.avg_out
                        new_degs[i,j,3] = stat.max_out
                    else
                        new_degs[i,j,1] = stat.min_in
                        new_degs[i,j,2] = stat.avg_in
                        new_degs[i,j,3] = stat.max_in
                    end
                end
            end

            data_label_in_colors = zeros(Float32, num_colors)
            for i in 1:num_colors
                data_label_in_colors[i] = 0
                for data_label in new_data_labels
                    if data_label == -1
                        data_label_in_colors[i] = 1
                        break
                    end
                    if data_label in summary.color_filters[new_color]
                        data_label_in_colors[i] = 1/summary.color_label_cardinality[j][new_label]
                    end
                end
            end
            old_color::Int16 = 1
            @turbo for i in 1:num_old_paths
                old_color = partial_paths[parent_idx, i]
                for j in 1:num_colors
                    new_partial_weights[1, j, i] = partial_weights[1, i] * new_degs[old_color, j, 1] * data_label_in_colors[j]
                    new_partial_weights[2, j, i] = partial_weights[2, i] * new_degs[old_color, j, 2] * data_label_in_colors[j]
                    new_partial_weights[3, j, i] = partial_weights[3, i] * new_degs[old_color, j, 3] * data_label_in_colors[j]
                end
            end

            nonzeros::Int = 0
            @turbo for i in 1:num_old_paths
                for j in 1:num_colors
                    nonzeros += new_partial_weights[2, j, i] > 0.0
                end
            end

            num_query_nodes = length(current_query_nodes)
            @turbo for i in 1:num_old_paths
                for j in 1:num_colors
                    for k in 1:num_query_nodes-1
                        new_partial_paths[k, j, i] = partial_paths[k, i]
                    end
                    new_partial_paths[num_query_nodes,  j, i] = j
                end
            end

            partial_paths = zeros(Color, num_query_nodes, nonzeros)
            partial_weights = zeros(Float64, 3, nonzeros)
            new_index = 1
            for i in 1:num_old_paths
                for j in 1:num_colors
                    @inbounds new_partial_weights[2, j, i] <= 0.0 && continue
                    for k in 1:num_query_nodes
                        @inbounds partial_paths[k, new_index]  = new_partial_paths[k, j, i]
                    end
                    @inbounds partial_weights[1, new_index]  = new_partial_weights[1, j, i]
                    @inbounds partial_weights[2, new_index]  = new_partial_weights[2, j, i]
                    @inbounds partial_weights[3, new_index]  = new_partial_weights[3, j, i]
                    new_index += 1
                end
            end
        else
            new_partial_paths = zeros(Color, length(current_query_nodes),  num_old_paths * num_colors)
            new_partial_weights = zeros(Float64, 3, num_old_paths * num_colors)

            # Update the partial paths using the parent-child combo that comes next from the query.
            for i in 1:num_old_paths
                old_color = partial_paths[parent_idx, i]
                # Account for colors with no outgoing children.
                if haskey(edge_deg, old_color)
                    for (new_color, degree_stats) in edge_deg[old_color]
                        # revamp the logic to use a set of labels rather than just one
                        # check if the data label(s) are in the color
                        data_label_in_color = false
                        for data_label in new_data_labels
                            if data_label == -1
                                data_label_in_color = true
                                continue
                            end
                            if data_label in summary.color_filters[new_color]
                                data_label_in_color = true
                            end
                        end
                        if !data_label_in_color
                            continue
                        end
                        degree_stats::DegreeStats = edge_deg[old_color][new_color]
                        for j in 1:length(current_query_nodes)-1
                            new_partial_paths[j, new_path_idx] = partial_paths[j, i]
                        end
                        new_partial_paths[length(current_query_nodes), new_path_idx] = new_color
                        if out_edge
                            new_partial_weights[1, new_path_idx] = partial_weights[1, i]*degree_stats.min_out
                            new_partial_weights[2, new_path_idx] = partial_weights[2, i]*degree_stats.avg_out
                            new_partial_weights[3, new_path_idx] = partial_weights[3, i]*degree_stats.max_out
                        else
                            new_partial_weights[1, new_path_idx] = partial_weights[1, i]*degree_stats.min_in
                            new_partial_weights[2, new_path_idx] = partial_weights[2, i]*degree_stats.avg_in
                            new_partial_weights[3, new_path_idx] = partial_weights[3, i]*degree_stats.max_in
                        end
                        if !(length(new_data_labels) == 1 && new_data_labels[1] == -1)
                            # we have already confirmed that the data label is in the color, but if the data label isn't -1
                            # then we need to scale down the result since we only want to consider one of the many nodes in the new color
                            # we also need to set the minimum to 0 but keep the maximum the same
                            new_partial_weights[1, new_path_idx] = 0
                            new_partial_weights[2, new_path_idx] /= summary.color_label_cardinality[new_color][new_label]
                        end
                        new_path_idx += 1
                    end
                end
            end
            partial_paths = new_partial_paths[:, 1:new_path_idx-1]
            partial_weights = new_partial_weights[:, 1:new_path_idx-1]
        end

        if (max_partial_paths !== nothing) && (size(partial_paths)[2] > max_partial_paths)
            partial_paths, partial_weights = sample_paths(partial_paths, partial_weights, max_partial_paths, sampling_strategy)
        end

        if (include_cycles)
            handle_extra_edges!(query, summary, partial_paths, partial_weights, current_query_nodes, visited_query_edges, usingStoredStats, only_shortest_path_cycle)
        end
    end

    # Sum over the calculated partial paths to get the final bounds.
    final_bounds = sum(partial_weights, dims=2)
    return final_bounds
end
