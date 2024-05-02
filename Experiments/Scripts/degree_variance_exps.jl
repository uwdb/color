using QuasiStableColors
include("../Experiments.jl")

# The goal of this file is to demonstrate the differences in degree as the number of colors increase

datasets = [human, aids, yeast, dblp]
partitioner = QuasiStable
partitioning_schemes = [[(partitioner, 0)],
                        [(partitioner, 3)],
                        [(partitioner, 15)],
                        [(partitioner, 31)],
                        [(partitioner, 63)],
                        [(partitioner, 127)]]
experiment_params = Vector{ExperimentParams}()

datapoint_datasets = []
num_colors = []
degree_deviations = []
for dataset in datasets
    println("Loading dataset: ", string(dataset))
    # obtain a copy of the graph
    data_graph = load_dataset(dataset)
    println("Processing dataset: ", string(dataset))
    for current_scheme in partitioning_schemes
        # color the graph
        println("Coloring: ", string(current_scheme))
        node_color_mapping = color_graph(data_graph, ColorSummaryParams(partitioning_scheme=current_scheme))
        println("Calculating Variance: ", string(current_scheme))
        # iterate through all colors and their vertices to figure out the standard deviation
        color_nodes_mapping::Dict{Int, Vector{Int}} = Dict(color=>[] for color in range(1, current_scheme[1][2]))
        for (node, color) in node_color_mapping
            current_list = get(color_nodes_mapping, color, [])
            push!(current_list, node)
            color_nodes_mapping[color] = current_list
        end
        # go through each color in the mapping and figure out the standard degree_deviations
        current_std_devs = []
        # iterate through each color used to partition the graph
        for color in keys(color_nodes_mapping)
            # find all the nodes that belong to the color
            nodes_in_color = get(color_nodes_mapping, color, [])
            color_neighbor_counts = Dict(possible_color => [0 for _ in nodes_in_color] for possible_color in keys(color_nodes_mapping)) # maps each (neighbor) color to a list of counts
            # start a dictionary mapping each color to a list of number of neighbors of that color
            # for each node, add its neighbor colors to the dictionary of overall neighbors by color
            for (i, node) in enumerate(nodes_in_color)
                current_count_per_color = counter(Int)
                current_node_neighbors = neighbors(data_graph.graph, node)
                # count all of the neighbors
                for neighbor in current_node_neighbors
                    color = node_color_mapping[neighbor]
                    color_neighbor_counts[color][i] += 1
                end
                # print(string(length(neighbors(data_graph.graph, node))) * ",")
            end
            # for each child color, get the std dev, then find the mean across all the colors
            # check if empty or 1 node in color or else std returns NaN
            std_devs_by_child_color = [length(neighbor_counts) <= 1 ? 0 : maximum(neighbor_counts)-minimum(neighbor_counts) for neighbor_counts in values(color_neighbor_counts)]
            current_std_dev = mean([std_dev  for std_dev in std_devs_by_child_color])
            push!(current_std_devs, current_std_dev)
        end
        println("Standard deviations: ", current_std_devs)
        overall_std_dev = mean(current_std_devs)
        push!(datapoint_datasets, dataset)
        push!(num_colors, current_scheme[1][2])
        push!(degree_deviations, isempty(current_std_devs) ? 0 : overall_std_dev)
    end
end
println("Datasets: ", string(datapoint_datasets))
println("Num Colors: ", string(num_colors))
println("Degree Deviations: ", string(degree_deviations))

# at this point, we have processed everything.

# save the resulting lists... not a csv because it's just a list of data points
filename = "degree_variance_results.txt"
destination = "Experiments/Results/"
results_file = open(destination * filename, "w")

println(results_file, "Datasets: ")
println(results_file, string(datapoint_datasets))

println(results_file, "Number of Colors: ")
println(results_file, string(num_colors))

println(results_file, "Degree Deviations: ")
println(results_file, string(degree_deviations))

close(results_file)

log_deviations = [deviation == 0 ? 0 : log10(deviation) for deviation in degree_deviations]

ENV["GKSwstype"]="100"
# now graph a scatter plot with lines connecting data points from the same dataset
p = plot(num_colors, log_deviations, group = datapoint_datasets, legend = :topright, size=(600, 400), linewidth=4, left_margin = 10mm, guidefont=14,xtickfont=12,ytickfont=12,legendfont=10,)
xlabel!(p, "Number of Colors")
ylabel!(p, "Degree Range log\$_{10}\$")
# savefig("Experiments/Results/Figures/degree_deviations_$(partitioner).pdf")
savefig("Experiments/Results/Figures/fig_2.pdf")