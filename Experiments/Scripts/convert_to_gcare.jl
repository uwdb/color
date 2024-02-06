include("../Experiments.jl")
const graph_id = 1; # consider replacing with the hash of the file

# take in a "subgraph matching" file, convert it to a data graph, then convert the data graph into a "gcare" file
function convert_dataset(path)
    data_graph = load_dataset(path, subgraph_matching_data=true)
    graph_file = open("Experiments/ConvertedGraphs/file.txt", "w")
    println(graph_file, "t # $graph_id")
    for node in 1:nv(data_graph.graph)
        current_node = node - 1
        current_labels = data_graph.vertex_labels[node]
        vertex_output = "v $current_node"
        for label in current_labels
            vertex_output *= " " * string(label)
        end
        println(graph_file, vertex_output)
    end
    for edge in edges(data_graph.graph)
        start_v = src(edge) - 1
        end_v = dst(edge) - 1
        label = only(data_graph.edge_labels[(start_v+1, end_v+1)])
        println(graph_file, "e $start_v $end_v $label")
    end
    close(graph_file)
end

convert_dataset("/homes/gws/dmbs/Cardinality-with-Colors/dataset/hprd/hprd.graph")