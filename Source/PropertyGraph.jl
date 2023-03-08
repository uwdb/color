using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
# Create a property graph struct which has the following things: underlying graph,
# and dictionaries which map edges to label sets and vertices to label sets

# queries only have one label per edge/node, data graphs can have any number.
struct PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Int, Dict{Int, Vector{Int}}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Dict{Int, Vector{Int}} # vertex_labels[n] = { labels }

    PropertyGraph(num_vertices::Int) = new(DiGraph(num_vertices), Dict(), Dict())

    PropertyGraph(g::DiGraph) = new(g, 
                                    Dict(src(e) => Dict(x => Vector{Int64}() for x in outneighbors(g, src(e))) for e in edges(g)),
                                    Dict(x => Vector{Int64}() for x in range(1, nv(g))))
end

struct QueryGraph
    graph::DiGraph
    edge_labels::Dict{Int, Dict{Int, Vector{Int}}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Dict{Int, Vector{Int}} # vertex_labels[n] = { labels }
    vertex_id_labels::Dict{Int, Int} # vertex_id_labels[n] = label

    QueryGraph(num_vertices::Int) = new(DiGraph(num_vertices), Dict(), Dict(), 
                                            Dict(x => -1 for x in range(1, num_vertices)))
    QueryGraph(g::DiGraph) = new(g, 
                                    Dict(src(e) => Dict(x => Vector{Int64}() for x in outneighbors(g, src(e))) for e in edges(g)),
                                    Dict(x => Vector{Int64}() for x in range(1, nv(g))),
                                    Dict(x => -1 for x in range(1, nv(g))))
end


# By default, node data labels are just their id - 1 in the G-Care benchmark.
# Assumes that data label matches the order the node was loaded.
function get_data_label(g::PropertyGraph, node::Int)
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
function update_node_labels(g::PropertyGraph, node::Int, node_labels::Vector{Int})
    g.vertex_labels[node] = node_labels
end

function update_node_labels(g::QueryGraph, node::Int, node_labels::Vector{Int})
    g.vertex_labels[node] = node_labels
end

# Adds a new edge to the graph then adds the new label to the edge's list of labels
function add_labeled_edge!(g::PropertyGraph, edge::Tuple{Int, Int}, edge_label::Int)
    add_edge!(g.graph, edge)
    parent_node = edge[1]
    child_node = edge[2]
    if !haskey(g.edge_labels, parent_node)
        g.edge_labels[parent_node] = Dict()
    end
    if !haskey(g.edge_labels[parent_node], child_node)
        g.edge_labels[parent_node][child_node] = []
    end
    if edge_label in g.edge_labels[parent_node][child_node]
        return
    end
    push!(g.edge_labels[parent_node][child_node], edge_label)
end 

function add_labeled_edge!(g::QueryGraph, edge::Tuple{Int, Int}, edge_label::Int)
    add_edge!(g.graph, edge)
    parent_node = edge[1]
    child_node = edge[2]
    if !haskey(g.edge_labels, parent_node)
        g.edge_labels[parent_node] = Dict()
    end
    if !haskey(g.edge_labels[parent_node], child_node)
        g.edge_labels[parent_node][child_node] = []
    end
    if edge_label in g.edge_labels[parent_node][child_node]
        return
    end
    push!(g.edge_labels[parent_node][child_node], edge_label)
end 