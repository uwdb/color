using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs

abstract type PropertyGraph end

struct DataGraph <: PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Tuple{Int, Int}, Vector{Int}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Vector{Vector{Int}} # vertex_labels[n] = { labels }

    DataGraph(num_vertices::Int) = DataGraph(DiGraph(num_vertices))

    DataGraph(g::DiGraph) = new(g, 
                                    Dict((src(e), dst(e)) => Vector{Int64}() for e in edges(g)),
                                    [Vector{Int}() for v in 1:nv(g)])
end

struct QueryGraph <: PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Tuple{Int, Int}, Vector{Int}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Vector{Int} # vertex_labels[n] = label
    vertex_id_labels::Vector{Int} # vertex_id_labels[n] = label

    QueryGraph(num_vertices::Int) = QueryGraph(DiGraph(num_vertices))
    QueryGraph(g::DiGraph) = new(g, 
                                    Dict((src(e), dst(e)) => Vector{Int64}() for e in edges(g)),
                                    [-1 for v in 1:nv(g)],
                                    [-1 for v in 1:nv(g)])
end

# By default, node data labels are just their id - 1 in the G-Care benchmark.
# Assumes that data label matches the order the node was loaded.
function get_data_label(g::DataGraph, node::Int)
    return node - 1
end

# For query graphs specifically, data labels can be arbitrary and must be kept track of.
function get_data_label(g::QueryGraph, node::Int)
    return g.vertex_id_labels[node]
end

# Only Query graphs need to have their data labels updated
function change_node_id!(g::QueryGraph, node::Int, id_label::Int)
    g.vertex_id_labels[node] = id_label
end

# Replaces the node's labels
function update_node_labels!(g::DataGraph, node::Int, node_labels::Vector{Int})
    g.vertex_labels[node] = node_labels
end

# Replaces the node's labels
function update_node_labels!(g::QueryGraph, node::Int, node_label::Int)
    g.vertex_labels[node] = node_label
end

# Replaces the node's data labels
function update_data_labels!(g::QueryGraph, node::Int, data_label::Int)
    g.vertex_id_labels[node] = data_label
end

# Adds a new edge to the graph then adds the new label to the edge's list of labels
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