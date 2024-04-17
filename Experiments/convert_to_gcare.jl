include("Experiments.jl")
# TODO for running on alley: lubm80, yeast, dblp, eu2005, patents, youtube
const graph_id = 1; # consider replacing with the hash of the file

# take in a "subgraph matching" file, convert it to a data graph, then convert the data graph into a "gcare" file
function convert_dataset(path, filename, destination; query=false)
    data_graph = load_dataset(path, subgraph_matching_data=true)
    graph_file = open(destination * filename * ".txt", "w")
    println(graph_file, "t # $graph_id")
    for node in 1:nv(data_graph.graph)
        current_node = node - 1
        current_labels = data_graph.vertex_labels[node]
        vertex_output = "v $current_node"
        
        for label in current_labels
            vertex_output *= " " * string(label)
        end
        if (query)
            vertex_output *= " " * string(-1)
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

# Convert graphs into gcare format

datasets = ["dblp", "eu2005", "patents", "youtube"]
# datasets =["yeast", ]
for dataset in datasets
    # convert data graphs
    # convert_dataset("/homes/gws/dmbs/Cardinality-with-Colors/dataset/" * dataset * "/" * dataset * ".graph", dataset,"Experiments/ConvertedGraphs/")

    # convert query graphs
    query_directory = "../queryset/" * dataset
    query_paths = readdir(query_directory, join=true)
    println("Processing queries for ", dataset)
    for query in query_paths
        filename = split(basename(query), ".")[1]
        convert_dataset(query, filename, "ConvertedGraphs/queryset/" * dataset * "/", query=true)
    end
end