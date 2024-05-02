using Plots.PlotMeasures
using Graphs
using Random
include("../Experiments.jl")

# The goal of this file is to investigate how build time is affected by size of graph

graph_sizes = [200000, 400000, 600000, 800000, 1000000]
build_time_sums::Vector{Float16} = [0,0,0,0,0]
num_edge_labels = 50
num_vertex_labels = 50
max_labels_per_v = 5
num_trials = 20

# using a dummy graph to get rid of @elapsed startup
println("doing dummy run")
dummy_graph = SimpleDiGraph(100, 100)
dummy_data = DataGraph(dummy_graph)
dummy_build_time = @elapsed generate_color_summary(dummy_data)

for i in 1:num_trials
    println("Trial: ", i)
    for j in eachindex(graph_sizes)
        # create a random graph with the desired total number of edges and vertices
        graph_size = graph_sizes[j]
        current_graph = SimpleDiGraph(convert(Int, graph_size / 2), convert(Int, graph_size / 2))
        built_graph = DataGraph(current_graph)
        # Randomly assign variables
        random_vertex_labels = [[rand(-1:num_vertex_labels) for num_v in 1:rand(1:max_labels_per_v)] for v in 1:nv(current_graph)]
        random_edge_labels = Dict()
        for edge in edges(current_graph)
            random_edge_labels[(src(edge), dst(edge))] = [rand(-1:num_edge_labels)]
        end
        built_graph.edge_labels = random_edge_labels
        built_graph.vertex_labels = random_vertex_labels
        # time the creation of the color summary
        build_time = @elapsed generate_color_summary(built_graph)
        build_time_sums[j] += build_time
    end
end
# find the average
average_build_times = build_time_sums ./ num_trials

ENV["GKSwstype"]="100"
p = bar(graph_sizes, 
    average_build_times,
    # thickness_scaling=1.25,
    bottom_margin = 20px,
    top_margin = 20px,
    left_margin = 10mm,
    size=(600, 400),
    titlefont = (12, :black),
    tickfont = (12, :black),
    guidefont = (15, :black),
    legend = false)
xlabel!("Graph Size (V+E)")
ylabel!("Build Time (s)")
savefig(p, "Experiments/Results/Figures/fig_10")

