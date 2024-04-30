using Graphs
"""
The ColorSummary struct holds statistical information associated with the colored graph.
It keeps detailed information about the number of edges between colors of a particular color and which land in
a particular color. Note that `-1` is used to represent a "wildcard" label. These do not appear in the data graph,
but they do occur in the query graph.
"""
mutable struct ColorSummary{DS}
    # for outdegrees, c2 is the color of the outneighbor
    # for indegrees, c2 is the color of the inneighbor
    # v2 represents the label of the node in c1
    color_label_cardinality::Dict{Color, Dict{Int, Int}} # color_label_cardinality[c][v] = num_vertices
    edge_deg::Dict{Int, Dict{Int, Dict{Color, Dict{Color, DS}}}} # edge_deg[e][v2][c1][c2] = degreestat
    color_filters::Dict{Color, SmallCuckoo} # color_filters[c] = filter
    color_full::Set{Color} # Denotes if a color's filter is full
    cycle_probabilities::Dict{CyclePathAndColors, Float64} # cycle_probabilities[[c1, c2], path] = likelihood
    cycle_length_probabilities::Dict{Int, Float64} #cycle_probabilities[path_length] = likelihood
    max_cycle_size::Int
    total_edges::Int
    total_nodes::Int
    num_colors::Int
    total_added_edges::Int
end

"""
Generates a color summary for the given DataGraph using specific parameters about the coloring method.
# Arguments
- g::DataGraph - the DataGraph to color
- params::ColorSummaryParams - parameters describing how to color the graph, including sampling and coloring techniques
- verbose - enables printed messages for debugging
- timing_vec::Vector{Float64} - a Vector used to store statistics about the time spent building the summary
"""
function generate_color_summary(g::DataGraph, params::ColorSummaryParams=ColorSummaryParams(); verbose=0, timing_vec::Vector{Float64} = Float64[])
    DS = params.deg_stats_type
    if (verbose > 0)
        println("Started coloring")
    end
    coloring_time = time()
    color_filters::Dict{Color, SmallCuckoo} = Dict()
    color_label_cardinality::Dict{Color, Any} = Dict()
    color_hash::Dict{NodeId, Color} = color_graph(g, params)
    num_colors = maximum(values(color_hash); init = 0)
    color_sizes = [0 for _ in 1:num_colors]
    for c in values(color_hash)
        color_sizes[c] += 1
    end
    coloring_time = time() - coloring_time
    push!(timing_vec, coloring_time)

    if (verbose > 0)
        println("Started cycle counting")
    end
    cycle_counting_time = time()
    cycle_probabilities::Dict{CyclePathAndColors, Float64} = Dict()
    cycle_length_probabilities::Dict{Int, Float64} = Dict()
    cycle_probabilities, cycle_length_probabilities = join_table_cycle_likelihoods(g, color_hash, params.max_cycle_size, params.max_partial_paths)
    cycle_counting_time = time() - cycle_counting_time
    push!(timing_vec, cycle_counting_time)

    current_color = 1;
    if (verbose > 0)
        println("Started bloom filters")
    end
    bloom_filter_time = time()
    for color in eachindex(color_sizes)
        num_nodes = max(1, color_sizes[color])
        accepted_error = 0.00001
        cuckoo_params = constrain(SmallCuckoo, fpr=accepted_error, capacity=num_nodes)
        color_filters[current_color] = SmallCuckoo{cuckoo_params.F}(cuckoo_params.nfingerprints)
        current_color += 1
    end
    bloom_filter_time = time() - bloom_filter_time
    push!(timing_vec, bloom_filter_time)


    if (verbose > 0)
        println("Started cardinality counts")
    end
    cardinality_counts_time = time()
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

        # increment counter for all labels, including wildcards
        inc!(color_label_cardinality[color], -1)
        for label in g.vertex_labels[node]
            inc!(color_label_cardinality[color], label)
        end
    end
    cardinality_counts_time = time() - cardinality_counts_time
    push!(timing_vec, cardinality_counts_time)

    if (verbose > 0)
        println("Started tracking statistics")
    end
    edge_stats_time = time()
    # We keep separate degree statistics for in-degree and out-degree.
    color_to_color_edge_list::Dict{Int, Dict{Int, Dict{Color, Dict{Color, Vector{Tuple{NodeId, NodeId, Bool}}}}}} = Dict()
    for x in vertices(g.graph)
        c1 = color_hash[x]
        for y in outneighbors(g.graph,x)
            c2 = color_hash[y]
            edge_labels = []
            copy!(edge_labels, g.edge_labels[(x, y)])
            push!(edge_labels, -1)
            y_vertex_labels = []
            copy!(y_vertex_labels, g.vertex_labels[y])
            push!(y_vertex_labels, -1)
            x_vertex_labels = []
            copy!(x_vertex_labels, g.vertex_labels[x])
            push!(x_vertex_labels, -1)

            # Since each edge/vertex can have multiple labels associated with it,
            # we must count each edge/vertex label separately in our counter.
            # Additionally, we need to include a count for the implicit `-1` wildcard label.
            for edge_label in edge_labels
                if !haskey(color_to_color_edge_list, edge_label)
                    color_to_color_edge_list[edge_label] = Dict()
                end
                # We also need to consider both outgoing and incoming labels
                for vertex_label in y_vertex_labels
                    if !haskey(color_to_color_edge_list[edge_label], vertex_label)
                        color_to_color_edge_list[edge_label][vertex_label] = Dict()
                    end
                    if !haskey(color_to_color_edge_list[edge_label][vertex_label], c1)
                        color_to_color_edge_list[edge_label][vertex_label][c1] = Dict()
                    end
                    if !haskey(color_to_color_edge_list[edge_label][vertex_label][c1], c2)
                        color_to_color_edge_list[edge_label][vertex_label][c1][c2] = []
                    end
                    push!(color_to_color_edge_list[edge_label][vertex_label][c1][c2], (x, y, true))
                end

                for vertex_label in x_vertex_labels
                    if !haskey(color_to_color_edge_list[edge_label], vertex_label)
                        color_to_color_edge_list[edge_label][vertex_label] = Dict()
                    end
                    if !haskey(color_to_color_edge_list[edge_label][vertex_label], c2)
                        color_to_color_edge_list[edge_label][vertex_label][c2] = Dict()
                    end
                    if !haskey(color_to_color_edge_list[edge_label][vertex_label][c2], c1)
                        color_to_color_edge_list[edge_label][vertex_label][c2][c1] = []
                    end
                    push!(color_to_color_edge_list[edge_label][vertex_label][c2][c1], (y, x, false))
                end
            end
        end
    end

    edge_deg::Dict{Int, Dict{Int, Dict{Color, Dict{Color, DegreeStats}}}} = Dict()
    for edge_label in keys(color_to_color_edge_list)
        edge_deg[edge_label] = Dict()
        for vertex_label in keys(color_to_color_edge_list[edge_label])
            edge_deg[edge_label][vertex_label] = Dict()
            for c1 in keys(color_to_color_edge_list[edge_label][vertex_label])
                edge_deg[edge_label][vertex_label][c1] = Dict()
                for c2 in keys(color_to_color_edge_list[edge_label][vertex_label][c1])
                    edge_deg[edge_label][vertex_label][c1][c2] = DS(g, color_to_color_edge_list[edge_label][vertex_label][c1][c2], color_sizes[c1])
                end
            end
        end
    end

    if (verbose > 0)
        println("Finished tracking statistics")
    end
    edge_stats_time = time() - edge_stats_time
    push!(timing_vec, edge_stats_time)

    return ColorSummary{DS}(color_label_cardinality, edge_deg, color_filters, Set{Color}(),
                cycle_probabilities, cycle_length_probabilities, params.max_cycle_size,
                 ne(g.graph), nv(g.graph), num_colors, 0)
end

"""
Approximates the probability of the cycle existing by using the degree into the landing node
and the total number of nodes in the landing node.
# Arguments
- edge_label - the label of the edge closing the cycle
- child_label - the label of the ending vertex of the closing edge
- parent_color - the color of the starting vertex of the closing edge
- child_color - the color of the ending vertex of the closing edge
- summary::ColorSummary - the graph summary used to calculate this likelihood
"""
@inline function get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary::ColorSummary)
    if (summary.color_label_cardinality[child_color][child_label] == 0)
        println("issue with independent cycle likelihood")
    end
    return get_out_deg_estimate(summary.edge_deg[edge_label][child_label][parent_color][child_color])/summary.color_label_cardinality[child_color][child_label]
end

"""
Approximates the probability of the cycle existing by treating each edge in the graph as an
independent event and considering the closure of the cycle as one such event.
When colors are unavailable, this coarser probability is used.
# Arguments
- summary::ColorSummary - the graph summary used to calculate this likelihood
"""
@inline function get_independent_cycle_likelihood(summary::ColorSummary)
    return min(1.0, 5.0 * summary.total_edges / summary.total_nodes ^ 2)
end


"""
Uses a colored graph to calculate the odds of cycles closing between nodes of specific colors.
Returns a cycle_likelihoods Dict mapping a cycle description (path directionality/length and colors) to its likelihood of existing,
and a more general cycle_length_likelihoods Dict mapping a cycle size to its likelihood of existing.
# Arguments
- g::DataGraph - the data graph used to generate statistics
- color_hash - a mapping of each node to its assigned color
- cycle_size::Int - the maximum cycle size to calculate and store statistics for
- max_partial_paths - the maximum number of partial paths to use during statistics calculating, enabling sampling
"""
function join_table_cycle_likelihoods(g::DataGraph, color_hash, cycle_size::Int, max_partial_paths)
    cycle_length_likelihoods::Dict{Int, Float64} = Dict()
    cycle_likelihoods::Dict{CyclePathAndColors, Float64} = Dict()
    if (cycle_size < 2)
        return cycle_likelihoods, cycle_length_likelihoods
    end

    # define a specific_edge as [n1, n2, c1, c2, d]
    # we can make a table mapping these to their counts, or just have an array of specific_edges

    # first, take each edge from the original graph and convert them into their new more specific forms.
    # we need to include the reversed version of the edge in case a path between nodes doesn't
    # start with a forward edge.
    # detailed_edges::Set{Tuple{Int, Int, Int, Int, Vector{Bool}}} = Set()

    println("Initializing Detailed Edges Dict: ", time())
    # map start_node => vector of edges with that start node
    detailed_edges::Dict{Int, Vector{Tuple{Int, Int, Int, Int, Vector{Bool}}}} = Dict(i => [] for i in vertices(g.graph))
    for edge in edges(g.graph)
        # detailed edge = [n1, n2, c1, c2, [d]]
        detailed_edge = (src(edge), dst(edge), color_hash[src(edge)], color_hash[dst(edge)], [true])
        detailed_reverse_edge = (dst(edge), src(edge), color_hash[dst(edge)], color_hash[src(edge)], [false])
        push!(detailed_edges[detailed_edge[1]], detailed_edge)
        push!(detailed_edges[detailed_reverse_edge[1]], detailed_reverse_edge)
    end

    # create tables for each size of cycle/path
    stored_cycles::Dict{CyclePathAndColors, Float32} = Dict() # this stores summary info representing the path lengths we want to close
    stored_paths::Dict{CyclePathAndColors, Float32} = Dict() # this stores summary info representing the path lengths that actually closed
    # summary info = [c1, c2, [d]]

    println("Initializing Cycle Path Calculation: ", time())
    # initialize with size two data
    # start up the "current joins" vector
    updated_paths::Dict{Tuple{Int, Int, Int, Int, Vector{Bool}}, Float64} = Dict() # stores our progress as we repeatedly join
    vertex_sample = collect(keys(detailed_edges))
    if length(vertex_sample) > max_partial_paths
        vertex_sample = sample(collect(keys(detailed_edges)), max_partial_paths; replace=false)
    end
    for vertex in vertex_sample
        for edge in detailed_edges[vertex]
            summary_info = CyclePathAndColors(edge[5], (edge[3], edge[4]))
            updated_paths[edge] = 1.0
            stored_paths[summary_info] = get(stored_paths, summary_info, 0) + 1.0
            edge_dir = edge[5][1]
            is_closed = (edge[1], edge[2], color_hash[edge[1]], color_hash[edge[2]], [!edge_dir]) in detailed_edges[edge[1]]
            stored_cycles[summary_info] = get(stored_cycles, summary_info, 0) + (is_closed ? 1.0 : 0.0)
        end
    end

    println("Performing Cycle Path Calculation: ", time())
    # for each cycle size...
    for current_cycle_size in 2: cycle_size
        println("Performing Cycle Size (", current_cycle_size, ") Calculation: ", time())
        # new_paths stores the current detailed paths and their count
        new_paths::Dict{Tuple{Int, Int, Int, Int, Vector{Bool}}, Float64} = Dict()
        joined_path::Tuple{Int, Int, Int, Int, Vector{Bool}} = (0,0,0,0,[])
        # join/aggregate all of the current paths with their connected edges
        # runtime of size-n iteration through the loop is o(e^[n-1]), where e is size of edges
        if (length(updated_paths) > max_partial_paths)
            new_dict::Dict{Tuple{Int, Int, Int, Int, Vector{Bool}}, Float64} = Dict()
            current_paths = collect(keys(updated_paths))
            sampled_paths = sample(current_paths, max_partial_paths; replace=false)
            sampled_weight = 0.0
            total_weight = sum(values(updated_paths))
            for path in sampled_paths
                sampled_weight += updated_paths[path]
            end
            for path in sampled_paths
                new_dict[path] = sampled_weight == 0 ? 0 : updated_paths[path] / sampled_weight * total_weight
            end
            updated_paths = new_dict
        end
        random_numbers = rand(UInt32, 5)
        for condensed_path in keys(updated_paths)
            # first join/aggregate everything
            num_neighbors = length(detailed_edges[condensed_path[2]])
            sampled_edges = detailed_edges[condensed_path[2]][(random_numbers .% num_neighbors) .+ 1]
#            sampled_edges = detailed_edges[condensed_path[2]]
            for edge in sampled_edges
                # [n1, n2, c1, c2, [d]]
                joined_path = (condensed_path[1], edge[2], condensed_path[3], edge[4], cat(condensed_path[5], edge[5], dims=1))
                # this aggregation should improve runtime...
                new_paths[joined_path] = get(new_paths, joined_path, 0) + updated_paths[condensed_path]*length(detailed_edges[condensed_path[2]])/length(sampled_edges)
            end
        end
        # go through all of the extended paths and track which ones close in a cycle
        for joined_path in keys(new_paths) # o(e)
            summary_info = CyclePathAndColors(joined_path[5], (joined_path[3], joined_path[4]))
            stored_paths[summary_info] = get(stored_paths, summary_info, 0) + new_paths[joined_path]
            # if the cycle exists, store the cycle in the overall cycles
            if (joined_path[1], joined_path[2], joined_path[3], joined_path[4], [false]) in detailed_edges[joined_path[1]]
                stored_cycles[summary_info] = get(stored_cycles, summary_info, 0) + new_paths[joined_path]
            end
        end
        # now, only look at the paths that were able to be joined
        updated_paths = new_paths
    end

    # now go through and aggregate all duplicates
    MinPathWeight = 25
    default_cycle_weights = Dict()
    default_cycle_counts = Dict()
    cycle_length_weights = Dict(i => 0.0 for i in 1:cycle_size)
    cycle_length_counts = Dict(i => 0.0 for i in 1:cycle_size)
    for path in keys(stored_paths)
        # CyclePathAndColors => count
        if stored_paths[path] >= MinPathWeight
            cycle_likelihoods[path] = get(stored_cycles, path, 0) / stored_paths[path]
        end
        default_path = color_path_to_default(path)
        default_cycle_weights[default_path] = get(default_cycle_weights, default_path, 0) + get(stored_cycles, path, 0)
        default_cycle_counts[default_path] = get(default_cycle_counts, default_path, 0) + stored_paths[path]

        path_length = length(path.path)
        cycle_length_weights[path_length] = get(cycle_length_weights, path_length, 0) + get(stored_cycles, path, 0)
        cycle_length_counts[path_length] = get(cycle_length_counts, path_length, 0) + stored_paths[path]
    end
    for default_path in keys(default_cycle_weights)
        if default_cycle_counts[default_path] < MinPathWeight
            continue
        end
        cycle_likelihoods[default_path] = default_cycle_weights[default_path] / default_cycle_counts[default_path]
    end
    cycle_length_likelihoods = Dict(i => cycle_length_counts[i] == 0 ? 0 : cycle_length_weights[i] / cycle_length_counts[i] for i in 1:cycle_size)
    return cycle_likelihoods, cycle_length_likelihoods
end
