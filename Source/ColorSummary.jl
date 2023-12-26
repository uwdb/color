using Graphs
using Probably

# The ColorSummary struct holds statistical information associated with the colored graph.
# It keeps detailed information about the number of edges between colors of a particular color and which land in
# a particular color. Note that `-1` is used to represent a "wildcard" label. These do not appear in the data graph,
# but they do occur in the query graph.
mutable struct ColorSummary{DS}
    # for outdegrees, c2 is the color of the outneighbor
    # for indegrees, c2 is the color of the inneighbor
    # v2 represents the label of the node in c1
    color_label_cardinality::Dict{Color, Dict{Int, Int}} # color_label_cardinality[c][v] = num_vertices
    edge_deg::Dict{Int, Dict{Int, Dict{Color, Dict{Color, DS}}}} # edge_deg[e][v2][c1][c2] = degreestat
    color_filters::Dict{Color, SmallCuckoo} # color_filters[c] = filter
    cycle_probabilities::Dict{CyclePathAndColors, Float64} # cycle_probabilities[[c1, c2], path] = likelihood
    cycle_length_probabilities::Dict{Int, Float64} #cycle_probabilities[path_length] = likelihood
    max_cycle_size::Int
    total_edges::Int
    total_nodes::Int
    num_colors::Int
    total_added_edges::Int
    added_color::Bool
end


function generate_color_summary(g::DataGraph, params::ColorSummaryParams=ColorSummaryParams(); verbose=0, timing_vec::Vector{Float64} = Float64[], use_cycle_join_table=true)
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
    if use_cycle_join_table
        cycle_probabilities, cycle_length_probabilities = join_table_cycle_likelihoods(g, color_hash, params.max_cycle_size, params.max_partial_paths)
    else
        cycle_probabilities, cycle_length_probabilities = get_color_cycle_likelihoods(params.max_cycle_size,
                                                                                            g,
                                                                                            color_hash,
                                                                                            params.max_partial_paths)
    end
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

    # println("color cardinality: ", color_label_cardinality)
    # println("edge deg: ", edge_deg)
    # println("number of nodes: ", nv(g.graph))
    # println("number of edges: ", ne(g.graph))

    return ColorSummary{DS}(color_label_cardinality, edge_deg, color_filters,
                cycle_probabilities, cycle_length_probabilities, params.max_cycle_size,
                 ne(g.graph), nv(g.graph), num_colors, 0, false)
end

function color_hash_to_groups(color_hash, num_colors)
    # the color hash maps from a node to its assigned color
    node_groups::Vector{Vector{Int}} = [[] for color in 1:num_colors]
    for node in keys(color_hash)
        push!(node_groups[color_hash[node]], node)
    end
    return node_groups
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
    numEntries = numEntries*3 # To account for min, avg, and max stats
    numEntries += length(summary.cycle_probabilities)
    return
end

# when we generate the color summary, we need to have a method  to calculate the odds of a cycle

# approximates the probability of the cycle existing by using the degree into the landing node
# and the total number of nodes in the landing node
@inline function get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary::ColorSummary)
    if (summary.color_label_cardinality[child_color][child_label] == 0)
        println("issue with independent cycle likelihood")
    end
    return get_out_deg_estimate(summary.edge_deg[edge_label][child_label][parent_color][child_color])/summary.color_label_cardinality[child_color][child_label]
end

# When colors are unavailable, we use a coarser probability
@inline function get_independent_cycle_likelihood(summary::ColorSummary)
    return min(1.0, 5.0 * summary.total_edges / summary.total_nodes ^ 2)
end


function join_table_cycle_likelihoods(g::DataGraph, color_hash, cycle_size::Int, max_partial_paths)
    # define a specific_edge as [n1, n2, c1, c2, d]
    # we can make a table mapping these to their counts, or just have an array of specific_edges

    # first, take each edge from the original graph and convert them into their new more specific forms.
    # we need to include the reversed version of the edge in case a path between nodes doesn't
    # start with a forward edge.
    # detailed_edges::Set{Tuple{Int, Int, Int, Int, Vector{Bool}}} = Set()

    # map start_node => vector of edges with that start node
    detailed_edges::Dict{Int, Set{Tuple{Int, Int, Int, Int, Vector{Bool}}}} = Dict(i => Set() for i in vertices(g.graph))
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

    # initialize with size two data
    # start up the "current joins" vector
    updated_paths::Dict{Tuple{Int, Int, Int, Int, Vector{Bool}}, Float64} = Dict() # stores our progress as we repeatedly join
    for edge_set in values(detailed_edges)
        for edge in edge_set
            summary_info = CyclePathAndColors(edge[5], (edge[3], edge[4]))
            updated_paths[edge] = 1.0
            stored_paths[summary_info] = 1.0
            stored_cycles[summary_info] = (edge[1], edge[2], color_hash[edge[1]], color_hash[edge[2]], [false]) in detailed_edges[edge[1]] ?
                                        1.0 : 0.0
        end
    end

    # for each cycle size...
    for current_cycle_size in 3: cycle_size
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
        for condensed_path in keys(updated_paths)
            # first join/aggregate everything
            for edge in detailed_edges[condensed_path[2]]
                # [n1, n2, c1, c2, [d]]
                joined_path = (condensed_path[1], edge[2], condensed_path[3], edge[4], cat(condensed_path[5], edge[5], dims=1))
                # this aggregation should improve runtime...
                new_paths[joined_path] = get(new_paths, joined_path, updated_paths[condensed_path]-1) + 1
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
    cycle_likelihoods::Dict = Dict()
    cycle_length_weights = Dict(i => 0.0 for i in 2:cycle_size)
    cycle_length_counts = Dict(i => 0.0 for i in 2:cycle_size)
    for path in keys(stored_paths)
        # CyclePathAndColors => count
        cycle_likelihoods[path] = get(stored_cycles, path, 0) / stored_paths[path]
        path_length = length(path.path) + 1
        cycle_length_weights[path_length] = get(cycle_length_weights, path_length, 0) + cycle_likelihoods[path]
        cycle_length_counts[path_length] = get(cycle_length_counts, path_length, 0) + 1
    end
    cycle_length_likelihoods = Dict(i => cycle_length_counts[i] == 0 ? 0 : cycle_length_weights[i] / cycle_length_counts[i] for i in 2:cycle_size)
    return cycle_likelihoods, cycle_length_likelihoods
end

# this is the one where we also have the directionality of the path
# returns a mapping from start/end-colors => cycle-likelihood
function get_color_cycle_likelihoods(max_size::Int, data::DataGraph, color_hash, max_partial_paths, min_partial_paths=50)
    # we map the path that needs to be closed to its likelihood
    # of actually closing use type-aliases (path = Vector{Bool})
    cycle_length_likelihoods::Dict{Int, Float64} = Dict()
    cycle_likelihoods::Dict{CyclePathAndColors, Float64} = Dict()
    if (max_size < 2)
        return cycle_likelihoods
    end
    default_color_pair::StartEndColorPair = (-1, -1)
    for i in 2:max_size
        i_path_weight, i_cycle_weight = (0.0,0.0)
        paths = generate_graphs(i - 1, 0, (Vector{DiGraph})([DiGraph(i)]), false)
        for path in paths
            total_path_weight, total_cycle_weight = (0.0,0.0)
            bool_graph = convert_path_graph_to_bools(path.graph)
            default_cycle_description = CyclePathAndColors(bool_graph, default_color_pair)
            likelihoods = approximate_color_cycle_likelihood(path, data, color_hash, max_partial_paths) # output a dictionary of start-end color pairs -> likelihood
            for color_pair in keys(likelihoods)
                if likelihoods[color_pair][1] == 0 || likelihoods[color_pair][2] == 0
                    continue
                end
                # if the number of paths is less than the minimum, we want to remove it from the likelihoods so it uses default instead
                if likelihoods[color_pair][1] >= min_partial_paths
                    current_cycle_description = CyclePathAndColors(bool_graph, color_pair)
                    # likelihoods[c1, c2] = [num_paths, num_cycles]
                    cycle_likelihoods[current_cycle_description] = (likelihoods[color_pair][1] == 0) ?
                                                                0 : (likelihoods[color_pair][2]) / likelihoods[color_pair][1]
                end
                total_path_weight += likelihoods[color_pair][1]
                total_cycle_weight += likelihoods[color_pair][2]
                i_path_weight += likelihoods[color_pair][1]
                i_cycle_weight += likelihoods[color_pair][2]
            end
            if total_path_weight > min_partial_paths/2
               cycle_likelihoods[default_cycle_description] = (total_path_weight == 0) ? 0 : total_cycle_weight/total_path_weight
            end
        end
        cycle_length_likelihoods[i] = (i_path_weight == 0) ? 0 : i_cycle_weight/i_path_weight
    end
    return cycle_likelihoods, cycle_length_likelihoods
end

# takes a path and converts it into a list of bools, each bool representing
# whether or not the next edge is going forwards
# ex: a => b => c <= d would become {true, true, false}
function convert_path_graph_to_bools(graph::DiGraph)
    # we say true = edge is going forwards
    bool_representation::Vector{Bool} = []
    for vertex in 1:(nv(graph)-1)
        push!(bool_representation, vertex in inneighbors(graph, vertex + 1))
    end
    return bool_representation
end


# returns a mapping of start/end-colors => path-count/cycle-count
function approximate_color_cycle_likelihood(path::QueryGraph, data::DataGraph, color_hash, max_partial_paths)
    nodes_to_keep = [1, nv(path.graph)]
    # add parameter for query nodes that shouldn't be aggregated out (the starting/ending nodes)
    # will now output a vector of tuples where first thing is a path (should only have start/end nodes after aggregations), second is the weight of the path
    # iterate through list to figure out cycle closure likelihoods
    partial_paths = get_subgraph_counts(path, data, max_partial_paths=max_partial_paths, nodes_to_keep=nodes_to_keep)
    color_matches::Dict{StartEndColorPair, Vector{Float64}} = Dict() # color_matches[c1,c2] = [num_paths, num_cycles]
    for path_and_weight in partial_paths
        # there should be only two nodes left in the path that aren't aggregated out
        # the returned subgraphs are actual paths, not color matches
        path = path_and_weight[1]
        path_weight = path_and_weight[2]
        # check the colors
        current_start_node = path[1]
        current_end_node = path[2]
        current_colors::StartEndColorPair = (color_hash[current_start_node], color_hash[current_end_node])
        if !(haskey(color_matches, current_colors))
            color_matches[current_colors] = [0, 0]
        end
        # if there is a closing edge, then count the entire weight of the path for the cycles as well
        # The path has a predefined directionality so we want to only find the likelihood that the last
        # node in the path wraps around to the beginning node
        cycle_weight = (in(current_end_node, inneighbors(data.graph, current_start_node))) ? path_weight : 0
        color_matches[current_colors][1] = color_matches[current_colors][1] + path_weight
        color_matches[current_colors][2] = color_matches[current_colors][2] + cycle_weight
    end

    # return mapping of color pairs -> path/cycle likelihood
    return color_matches
end

# For a given number of edges, generates all possible directed graphs
# with the given number of edges
function generate_graphs(desiredEdges::Int, finishedEdges::Int, graphs::Vector{DiGraph}, isCyclic::Bool)
    if (finishedEdges >= desiredEdges)
        # now we are guaranteed to have all graphs have the correct number of edges
        # and close the loop, so we can return this final result...
        if (isCyclic)
            graphs = filter!(g->ne(g)!=1, graphs)
        end
        query_graphs::Vector{QueryGraph} = []
        for g in graphs
            current_query = QueryGraph(g)
            for edge in edges(g)
                update_edge_labels!(current_query, (src(edge), dst(edge)), [-1])
            end
            push!(query_graphs, current_query)
        end
        return query_graphs
    end
    startNode = finishedEdges + 1 # the next edge starts at the end of the finished edge
    nextNode::Int = startNode + 1
    if desiredEdges == (finishedEdges + 1)
        # if the graph is cyclic, the last edge should go back to the starting node
        # if the graph isn't cyclic, the last edge should go to the last node remaining
        nextNode = isCyclic ? 1 : nextNode
    end
    newGraphs::Vector{DiGraph} = []
    # for each graph, add one edge going to the next node
    # or one edge coming from the next node
    for graph in graphs
        graphWithForwardEdge = copy(graph)
        add_edge!(graphWithForwardEdge, startNode, nextNode)
        push!(newGraphs, graphWithForwardEdge)
        # if we're on the last edge of a cycle, don't do this since it results in duplicate graphs
        if (!(isCyclic && desiredEdges == finishedEdges + 1))
            graphWithBackEdge = copy(graph)
            add_edge!(graphWithBackEdge, nextNode, startNode)
            push!(newGraphs, graphWithBackEdge)
        end
    end
    return generate_graphs(desiredEdges, finishedEdges + 1, newGraphs, isCyclic)
end
