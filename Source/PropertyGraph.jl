using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using Graphs
# Create a property graph struct which has the following things: underlying graph,
# and dictionaries which map edges to label sets and vertices to label sets

# queries only have one label per edge/node, data graphs can have any number.
struct PropertyGraph
    graph::DiGraph
    edge_labels::Dict{Int, Dict{Int, Array{Int}}} # edge_labels[n1][n2] = { labels }
    vertex_labels::Dict{Int, Array{Int}} # vertex_labels[n] = { labels }

    PropertyGraph(num_vertices::Int) = new(DiGraph(num_vertices), Dict(), Dict())
end


function add_labeled_node!(g::PropertyGraph, node::Int, node_labels::Array{Int})
    g.vertex_labels[node] = node_labels
end

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