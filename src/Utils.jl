"""
Creates a spanning tree and replaces the graph in the given Query.
# Arguments
- query::QueryGraph - the graph to find a spanning tree for.
"""
function get_spanning_tree(query::QueryGraph)
    # can keep the original QueryGraph object to preserve the labels, we just need
    # to change out query.graph with a spanning tree version.
    current_graph = query.graph
    undirected_graph = SimpleGraph(current_graph)
    undirected_spanning_tree = bfs_tree(undirected_graph, 1)  # Start a bfs tree from node 1
    # for acyclic queries, this will end up being the same graph since all the spanning trees in an
    # acyclic graph are the same graph
    new_graph = SimpleDiGraph(nv(query.graph))
    for edge in edges(undirected_spanning_tree)
        if has_edge(current_graph, src(edge), dst(edge))
            add_edge!(new_graph, src(edge), dst(edge))
        else
            add_edge!(new_graph, dst(edge), src(edge))
        end
    end
    query.graph = new_graph
end

"""
Finds the node order of the given graph that gives the minimum tree width and returns it as a list of nodes (also optionally returns the width).
# Arguments
- g::DiGraph - the graph to find a node order for.
- max_exact_width - the maximum exact width allowed for the node order.
- return_width - whether or not to also return the tree width of the resulting node order.
"""
function get_min_width_node_order(g::DiGraph; max_exact_width = 1, return_width=false)
    if (nv(g) < max_exact_width)
        nodes_processed = 1
        partial_orders::Dict{Set{Int32}, Tuple{Vector{Int}, Int32}} = Dict(Set([x]) =>([x], 0) for x in vertices(g))
        while nodes_processed < nv(g)
            for node_set in keys(partial_orders)
                best_order, best_width = partial_orders[node_set]
                neighbor_nodes = Set{Int32}()
                for existing_node in node_set
                    for neighbor in all_neighbors(g, existing_node)
                        if !(neighbor in node_set)
                            push!(neighbor_nodes, neighbor)
                        end
                    end
                end
                for next_node in neighbor_nodes
                    new_node_set::Set{Int32} = Set([next_node, node_set...])
                    new_order::Vector{Int32} = [best_order..., next_node]
                    new_width = 0
                    for v in new_order
                        if ! all([x in new_order for x in all_neighbors(g, v)])
                            new_width += 1
                        end
                    end
                    new_width = max(best_width, new_width)
                    if haskey(partial_orders, new_node_set)
                        if partial_orders[new_node_set][2] > new_width
                            partial_orders[new_node_set] = (new_order, new_width)
                        end
                    else
                        partial_orders[new_node_set] = (new_order, new_width)
                    end
                end
            end
            nodes_processed += 1
            for node_set in keys(partial_orders)
                if length(node_set) < nodes_processed
                    delete!(partial_orders, node_set)
                end
            end
        end
        min_width = minimum([x[2] for x in values(partial_orders)])
        for node_order_and_width in values(partial_orders)
            if node_order_and_width[2] == min_width
                return_width && return node_order_and_width[2]
                return node_order_and_width[1]
            end
        end
    else #NOTE THIS IS A PERFORMANCE BOTTLENECK
        min_width = nv(g)
        min_order::Vector{Int} = []
        for starting_node in vertices(g)
            max_width = 0
            visited_nodes::Vector{Int}= [starting_node]
            neighbor_count = Dict{Int, Int}(starting_node => degree(g, starting_node))
            living_neighbors = Set{Int}(all_neighbors(g, starting_node))
            current_width = 1
            while length(visited_nodes) < nv(g)
                new_width = nv(g)
                next_node = -1
                for potential_node in living_neighbors
                    new_neighbors = 0
                    killed_neighbors = 0
                    for potential_new_neighbor in all_neighbors(g, potential_node)
                        if !(potential_new_neighbor in visited_nodes)
                            new_neighbors = 1
                        elseif neighbor_count[potential_new_neighbor] == 1
                            killed_neighbors += 1
                        end
                    end
                    if current_width + new_neighbors - killed_neighbors <= new_width
                        next_node = potential_node
                        new_width = current_width + new_neighbors - killed_neighbors
                    end
                end

                neighbor_count[next_node] = 0
                for neighbor in all_neighbors(g, next_node)
                    if neighbor in visited_nodes
                        neighbor_count[neighbor] -= 1
                    else
                        neighbor_count[next_node] += 1
                        push!(living_neighbors, neighbor)
                    end
                end
                delete!(living_neighbors, next_node)
                push!(visited_nodes, next_node)
                current_width = new_width
                max_width = max(max_width, new_width)
            end
            if max_width <= min_width
                min_order = visited_nodes
                min_width = max_width
            end
        end
        return_width && return min_width
        return min_order
    end
end
