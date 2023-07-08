include("PropertyGraph.jl")
using AutoHashEquals

BoolPath = Vector{Bool}
StartEndColorPair = Vector{Int}
abstract type Comparable end
import Base .==
function ==(a::T, b::T) where T <: Comparable
    (a.path == b.path) && (a.colors == b.colors)
end
@auto_hash_equals mutable struct CyclePathAndColors
    path::BoolPath
    colors::StartEndColorPair
end

# The ColorSummary struct holds statistical information associated with the colored graph.
# It keeps detailed information about the number of edges between colors of a particular color and which land in
# a particular color. Note that `-1` is used to represent a "wildcard" label. These do not appear in the data graph,
# but they do occur in the query graph.
struct ColorSummary
    color_label_cardinality::Dict{Int, Dict{Int, Int}} # color_label_cardinality[c][v] = num_vertices
    edge_min_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_min_out_deg[e][v2][c1][c2] = min
    edge_min_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_min_in_deg[e][v2][c1][c2] = min
    edge_avg_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_avg_out_deg[e][v2][c1][c2] = avg
    edge_avg_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_avg_in_deg[e][v2][c1][c2] = avg
    edge_max_out_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_max_out_deg[e][v2][c1][c2] = max
    edge_max_in_deg::Dict{Int, Dict{Int, Dict{Int, Dict{Int, Float64}}}} # edge_max_in_deg[e][v2][c1][c2] = max
    color_filters::Dict{Int, BloomFilter} # color_filters[c] = filter
    cycle_probabilities::Dict{CyclePathAndColors, Float64} # cycle_probabilities[[c1, c2], path] = likelihood
    total_edges::Int
    total_nodes::Int
    # cycle_probabilities::Dict{Int, Float64}
    # for outdegrees, c2 is the color of the outneighbor
    # for indegrees, c2 is the color of the inneighbor
    # v2 represents the label of the node in c1
end

function generate_color_summary(g::DataGraph, numColors::Int; weighting=true, verbose=false, 
                                max_size=4, num_sample_nodes=1000, partitioner ="QuasiStable")
    color_hash = nothing
    color_sizes = [0 for _ in 1:numColors]
    color_filters = Dict()
    color_cardinality = Dict()
    color_label_cardinality = Dict()
    if partitioner == "QuasiStable"
        QSC = QuasiStableColors
        if (verbose) 
            println("Started coloring")
        end
        C = QSC.q_color(g.graph, n_colors=numColors, weighting=weighting)
        if (verbose)
            println("Finished coloring")
        end
        color_hash = QSC.node_map(C)
        color_sizes = [only(size(C.partitions[i])) for i in 1:length(C.partitions)] # EXTREMELY SUS
    elseif partitioner == "Hash"
        color_hash = Dict()
        for i in 1:nv(g.graph)
            color_hash[i] = (hash(i) % numColors) + 1
            color_sizes[color_hash[i]] += 1
        end
    elseif partitioner == "Degree"
        color_hash = Dict()
        degrees = sort(degree(g.graph))
        bucket_right_edges = []
        for i in 1:numColors
            degree_quantile = degrees[Int(floor(float(i)/numColors *length(degrees)))]
            if i > 1 && bucket_right_edges[i-1] >= degree_quantile
                push!(bucket_right_edges, bucket_right_edges[i-1] + 1)
            else
                push!(bucket_right_edges, degree_quantile)
            end
        end
        for i in 1:nv(g.graph)
            node_degree = degree(g.graph, i)
            for j in 1:numColors
                if node_degree <= bucket_right_edges[j]
                    color_hash[i] = j
                    color_sizes[color_hash[i]] += 1
                    break
                end
            end
        end

    elseif partitioner == "DirectedDegree"
        color_hash = Dict()
        num_degree_buckets = Int(floor(float(numColors)^.5))
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
            color_sizes[color_hash[i]] += 1
        end
    end
    
    # cycle_probabilities::Dict{Int, Dict{Vector{Bool}, Float64}} = get_cycle_likelihoods(max_size, g, num_sample_nodes)
    cycle_probabilities::Dict{CyclePathAndColors, Float64} = get_color_cycle_likelihoods(max_size, g, color_hash, num_sample_nodes)
    # the color hash should map node => color
    # have a paramter for the color hash, make changes to get_cycle_likelihoods

    # initialize color filters for data labels
    current_color = 1;
    if (verbose)
        println("Started bloom filters")
    end
    for color in eachindex(color_sizes)
        num_nodes = max(1, color_sizes[color])
        accepted_error = 0.00001
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
                                    edge_avg_out_deg, edge_avg_in_deg, edge_max_out_deg, edge_max_in_deg, color_filters,
                                     cycle_probabilities, ne(g.graph), nv(g.graph))
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
function get_independent_cycle_likelihood(edge_label, child_label, parent_color, child_color, summary::ColorSummary)
    summary.edge_avg_out_deg[edge_label][child_label][parent_color][child_color]/summary.color_label_cardinality[child_color][child_label]
end

# approximates the probability of a cycle existing based on the starting color of the path to be closed and the
# directionality of the path that will be closed
function get_start_color_cycle_likelihoods(max_cycle_size::Int, data::DataGraph, color_hash; num_samples_per_color::Int=0)
    # right now the color hash is a mapping of node -> color, but we want to invert that:
    color_nodes_mapping::Dict{Int, Vector{Int}} = Dict()
    for node in keys(color_hash)
        if (!haskey(color_nodes_mapping, color_hash[node]))
            color_nodes_mapping[color_hash[node]] = []
        end
        push!(color_nodes_mapping[color_hash[node]], node)
    end

    cycle_likelihoods::Dict{Int, Dict{Vector{Bool}, Float64}} = Dict() # [c1][bool_path] = likelihood
    # basically the same as the other method, except the query graph should now have a set of data labels attached to it

    for color in keys(color_nodes_mapping)
        cycle_likelihoods[color] = Dict()
        for size in 2: max_cycle_size
            paths = generate_graphs(size - 1, 0, (Vector{DiGraph})([DiGraph(size)]), false)
            for path_query in paths
                cycle = copy(path_query.graph)
                add_edge!(cycle, nv(cycle), 1)
                if (ne(cycle) == 1)
                    # accounts for case where the two edges from one node
                    # point to the same destination node
                    continue
                end
                cycle_query = QueryGraph(cycle)
                current_starting_nodes = nothing
                if !(num_samples_per_color === 0)
                    sample_size = min(num_samples_per_color, length(color_nodes_mapping[color]))
                    current_starting_nodes = sample(color_nodes_mapping[color], sample_size, replace=false)
                end
                numCycles::Float64 = get_exact_size(cycle_query, data, starting_nodes=current_starting_nodes)
                numPaths::Float64 = get_exact_size(path_query, data, starting_nodes=current_starting_nodes)
                bool_representation = convert_path_graph_to_bools(path_query.graph)
                likelihood = numCycles / numPaths
                cycle_likelihoods[color][bool_representation] = likelihood
            end
        end
    end
    return cycle_likelihoods
end

# gets the directed path from the start to finish node
function get_matching_graph(start::Int, finish::Int, query::QueryGraph)
    # remove the edge closing the cycle
    directed_graph_copy = copy(query.graph)
    rem_edge!(directed_graph_copy, start, finish)

    # convert the graph to be undirected
    graph_copy = Graph(copy(directed_graph_copy))
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

# approximates the probability of a cycle existing based on the directionality of the path that will be closed
function get_cycle_likelihoods(max_size::Int, data::DataGraph, num_sample_nodes)
    # we map the path that needs to be closed to its likelihood
    # of actually closing
    cycle_likelihoods::Dict{Vector{Bool}, Float64} = Dict()
    if (max_size < 2)
        return cycle_likelihoods
    end
    for i in 2:max_size
        println("Generating Cycles of Size: ", i)
        paths = generate_graphs(i - 1, 0, (Vector{DiGraph})([DiGraph(i)]), false)
        for path in paths
            cycle = copy(path.graph)
            add_edge!(cycle, nv(cycle), 1)
            if (ne(cycle) == 1)
                # accounts for case where the two edges from one node
                # point to the same destination node
                continue
            end
            cycleGraph = QueryGraph(cycle)
            likelihood = approximate_cycle_likelihood(path, cycleGraph, data, num_sample_nodes) # output a dictionary of start-end color pairs -> likelihood
            cycle_likelihoods[convert_path_graph_to_bools(path.graph)] = likelihood
        end
    end
    return cycle_likelihoods
end

# returns a mapping from start/end-colors => cycle-likelihood
function get_color_cycle_likelihoods(max_size::Int, data::DataGraph, color_hash, num_sample_nodes)
    # we map the path that needs to be closed to its likelihood
    # of actually closing
    # use type-aliases (path = Vector{Bool})
    cycle_likelihoods::Dict{CyclePathAndColors, Float64} = Dict()
    if (max_size < 2)
        return cycle_likelihoods
    end
    for i in 2:max_size
        paths = generate_graphs(i - 1, 0, (Vector{DiGraph})([DiGraph(i)]), false)
        for path in paths
            likelihoods = approximate_color_cycle_likelihood(path, data, color_hash, num_sample_nodes) # output a dictionary of start-end color pairs -> likelihood
            for color_pair in keys(likelihoods)
                bool_graph = convert_path_graph_to_bools(path.graph)
                current_cycle_description = CyclePathAndColors(bool_graph, color_pair)
                # likelihoods[c1, c2] = [num_paths, num_cycles]
                cycle_likelihoods[current_cycle_description] = (likelihoods[color_pair][1] == 0) ?
                                                            0 : (likelihoods[color_pair][2]) / likelihoods[color_pair][1]
            end
        end
    end
    return cycle_likelihoods
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

# for a specific path, calculates the
# probability that the cycle closes
function approximate_cycle_likelihood(path::QueryGraph, cycle::QueryGraph, data::DataGraph, num_sample_nodes)
    sampled_starting_nodes = nothing
    if !(num_sample_nodes === nothing)
        sampled_starting_nodes = sample(1:nv(data.graph), num_sample_nodes, replace=false)
    end
    numCycles::Float64 = get_exact_size(cycle, data, starting_nodes=sampled_starting_nodes) 
    # add parameter for query nodes that shouldn't be aggregated out (the starting/ending nodes)
    # will now output a vector of tuples where first thing is a path (should only have start/end nodes after aggregations), second is the weight of the path
    # iterate through list to figure out cycle closure likelihoods
    numPaths::Float64 = get_exact_size(path, data, starting_nodes=sampled_starting_nodes)
    # only find paths, use data graph to find closing edge if existing
    return numPaths != 0 ? numCycles / numPaths : 0
end

# returns a mapping of start/end-colors => path-count/cycle-count
function approximate_color_cycle_likelihood(path::QueryGraph, data::DataGraph, color_hash, num_sample_nodes)
    sampled_starting_nodes = nothing
    if !(num_sample_nodes === nothing)
        num_samples = min(nv(data.graph), num_sample_nodes)
        sampled_starting_nodes = sample(1:nv(data.graph), num_samples, replace=false)
    end
    nodes_to_keep = [1, nv(path.graph)]
    # add parameter for query nodes that shouldn't be aggregated out (the starting/ending nodes)
    # will now output a vector of tuples where first thing is a path (should only have start/end nodes after aggregations), second is the weight of the path
    # iterate through list to figure out cycle closure likelihoods
    partial_paths = get_subgraph_counts(path, data, starting_nodes=sampled_starting_nodes, nodes_to_keep=nodes_to_keep)
    color_matches::Dict{StartEndColorPair, Vector{Float64}} = Dict() # color_matches[c1,c2] = [num_paths, num_cycles]
    for path_and_weight in partial_paths
        # there should be only two nodes left in the path that aren't aggregated out
        # the returned subgraphs are actual paths, not color matches
        path = path_and_weight[1]
        path_weight = path_and_weight[2]
        # check the colors
        current_start_node = path[1]
        current_end_node = path[2]
        current_colors::StartEndColorPair = [color_hash[current_start_node], color_hash[current_end_node]]
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

# Given a cycle size, find the probability that a chain with the same number of nodes
# will also have an edge closing the cycle
# Note: used for old version that mapped cycle-size -> cycle-likelihood
function approximate_cycle_likelihood(max_cycle_size::Int, data::DataGraph)
    cycle_likelihoods::Dict{Int,Float64}
    for size in 2:max_cycle_size
        # convert to undirected graphs so we don't have to deal with duplicate generated graphs
        undirected_data = DataGraph(DiGraph(Graph(data.graph)))
        undirected_path = QueryGraph(DiGraph(path_graph(max_cycle_size)))
        undirected_cycle = QueryGraph(DiGraph(cycle_graph(max_cycle_size)))
        num_paths = get_exact_size(undirected_path, undirected_data)
        num_cycles = get_exact_size(undirected_cycle, undirected_data)
        cycle_likelihoods[size] = num_paths == 0 ? 0 : num_cycles / num_paths
    end
    return cycle_likelihoods
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

# Returns the shortest distance between two nodes not including their directly-connecting edge
# Note: previously used for path-length -> cycle-likelihood mapping
function shortestPathNotDirectlyConnected(startNode::Int, endNode::Int, query::QueryGraph)
    copiedDiGraph = copy(query.graph)
    # remove the edge we are trying to close in the query
    rem_edge!(copiedDiGraph, startNode, endNode)
    # find the shortest path between the start and end node
    # start by removing directionality because we're just concerned with the shortest
    # number of connections between the two nodes
    undirectedGraph = Graph(copiedDiGraph)
    ds = dijkstra_shortest_paths(undirectedGraph, startNode)
    return ds.dists[endNode]
end