
abstract type PropertyGraph end

"""
Representation of a Data Graph.
Data Graphs are unique because each node can have multiple labels but they have a unique assigned data labels.
"""
mutable struct DataGraph <: PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Tuple{Int, Int}, Vector{Int}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Vector{Vector{Int}} # vertex_labels[n] = { labels }

    DataGraph(num_vertices::Int) = DataGraph(DiGraph(num_vertices))

    DataGraph(g::DiGraph) = new(g,
                                    Dict((src(e), dst(e)) => Vector{Int64}([-1]) for e in edges(g)),
                                    [Vector{Int}([-1]) for v in 1:nv(g)])
end

"""
Representation of a Query Graph.
Query Graphs are unique because each node can only have one label and their data labels are arbitrary.
"""
mutable struct QueryGraph <: PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Tuple{Int, Int}, Vector{Int}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Vector{Vector{Int}} # vertex_labels[n] = labels
    vertex_id_labels::Vector{Vector{Int}} # vertex_id_labels[n] = labels

    QueryGraph(num_vertices::Int) = QueryGraph(DiGraph(num_vertices))
    QueryGraph(g::DiGraph) = new(g,
                                    Dict((src(e), dst(e)) => Vector{Int64}([-1]) for e in edges(g)),
                                    [[-1] for v in 1:nv(g)],
                                    [[-1] for v in 1:nv(g)])
end

"""
Returns the data label of the node in the Data Graph. Note that by default,
node data labels are just their ID - 1 in the G-Care benchmark.
Assumes that the data label matches the order that the node was loaded.
# Arguments
- g::DataGraph - the Data Graph that the node belongs to.
- node::Int - the vertex ID of the node that we find the data label for.
"""
function get_data_label(g::DataGraph, node::Int)
    return node - 1
end

"""
Obtains the data label of the node in the Query Graph.
For query graphs specifically, data labels can be arbitrary and must be tracked.
# Arguments
- g::QueryGraph - the Query Graph that the node belongs to.
- node::Int - the vertex ID of the node that we find the data label for.
"""
function get_data_label(g::QueryGraph, node::Int)
    return g.vertex_id_labels[node]
end

"""
Changes the data label of the node in the Query Graph.
Only Query Graphs have arbitrary data labels for specific nodes that
may need to be changed.
# Arguments
- g::QueryGraph - the Query Graph that the node belongs to.
- node::Int - the vertex ID of the node that we update the data label for.
- id_label::Int - the new data label to assign to the node.
"""
function change_node_id!(g::QueryGraph, node::Int, id_label::Int)
    g.vertex_id_labels[node] = [id_label]
end

"""
Replaces all of the given node's labels in the Data Graph.
For Data Graphs, nodes can have any number of labels.
# Arguments
- g::DataGraph - the Data Graph that the node belongs to.
- node::Int - the vertex ID of the node that we update the label for.
- node_labels::Vector{Int} - the new labels to use for the node.
"""
function update_node_labels!(g::DataGraph, node::Int, node_labels::Vector{Int})
    g.vertex_labels[node] = node_labels
end

"""
Replaces the given node's label in the Query Graph.
For Query Graphs, nodes only have one label.
# Arguments
- g::QueryGraph - the Query Graph that the node belongs to.
- node::Int - the vertex ID of the node that we update the label for.
- node_label::Int - the new label to use for the node.
"""
function update_node_labels!(g::QueryGraph, node::Int, node_label::Int)
    g.vertex_labels[node] = [node_label]
end

"""
Replaces the given node's data label in the Query Graph with one label.
# Arguments
- g::QueryGraph - the Query Graph that the node belongs to.
- node::Int - the vertex ID of the node that we update the label for.
- data_label::Int - the new data label to use for the node.
"""
function update_data_labels!(g::QueryGraph, node::Int, data_label::Int)
    g.vertex_id_labels[node] = [data_label]
end

"""
Replaces the given node's data label in the Query Graph with multiple label.
Used in cases where one query has different options for a specific node.
# Arguments
- g::QueryGraph - the Query Graph that the node belongs to.
- node::Int - the vertex ID of the node that we update the label for.
- data_labels::Vector{Int} - the new data labels to use for the node.
"""
function update_data_labels!(g::QueryGraph, node::Int, data_labels::Vector{Int})
    g.vertex_id_labels[node] = data_labels
end

"""
Adds a new edge to the given graph then adds the new label to the edge's list of labels.
If an edge between the vertices already exists, no updates will be made.
# Arguments
- g::PropertyGraph - the graph to update.
- edge::Tuple{Int, Int} - the edge connection (start, end) to add to the graph.
- edge_label::Int - the label to assign to the edge.
"""
function add_labeled_edge!(g::PropertyGraph, edge::Tuple{Int, Int}, edge_label::Int)
    add_edge!(g.graph, edge)
    if !haskey(g.edge_labels, edge)
        g.edge_labels[edge] = Vector{Int64}([])
    end
    if edge_label in g.edge_labels[edge]
        return
    end
    push!(g.edge_labels[edge], edge_label)
end

"""
Replaces the edge label(s) in the Query Graph for a specific edge.
# Arguments
- g::QueryGraph - the Query Graph to update
- edge::Tuple{Int, Int} - the edge connection (start, end) where the label(s) should be updated.
- edge_labels::Vector{Int} - a list of labels to assign to the edge.
"""
function update_edge_labels!(g::QueryGraph, edge::Tuple{Int, Int}, edge_labels::Vector{Int})
    g.edge_labels[edge] = edge_labels
end

"""
Adds a new (labeled) node to the Data Graph
# Arguments
- g::DataGraph - the Data Graph to update_data_labels
- node_labels - the labels to assign to the new node, if any.
"""
function add_labeled_node!(g::DataGraph, node_labels=[])
    add_vertex!(g.graph)
    push!(g.vertex_labels, node_labels)
end
