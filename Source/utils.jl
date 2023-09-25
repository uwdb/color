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


function get_min_width_node_order(g::DiGraph)
    if (nv(g) < 10)
        nodes_processed = 1
        partial_orders::Dict{Set{Int64}, Tuple{Vector{Int64}, Int64}} = Dict(Set([x]) =>([x], 0) for x in vertices(g))
        while nodes_processed < nv(g)
            for node_set in keys(partial_orders)
                best_order, best_width = partial_orders[node_set]
                neighbor_nodes = Set()
                for existing_node in node_set
                    for neighbor in all_neighbors(g, existing_node)
                        if !(neighbor in node_set)
                            push!(neighbor_nodes, neighbor)
                        end
                    end
                end
                for next_node in neighbor_nodes
                    new_node_set::Set{Int32} = Set([next_node])
                    new_node_set = union(new_node_set, node_set)
                    new_order = []
                    copy!(new_order, best_order)
                    push!(new_order, next_node)
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
                return node_order_and_width[1]
            end
        end
    else
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
end
