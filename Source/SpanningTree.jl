include("QuasiStableCardinalityEstimator.jl")

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