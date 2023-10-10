using Plots.PlotMeasures 
include("Experiments/build_color_summaries.jl")
include("Experiments/get_true_cardinalities.jl")
include("Experiments/load_datasets.jl")
include("Experiments/load_querysets.jl")
include("Experiments/run_estimators.jl")
include("Experiments/graph_results.jl")
include("Experiments/utils.jl")

# QuasiStable Hash Degree DirectedDegree SimpleLabel InOut LabelInOut NeighborEdges MostNeighbors

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=aids, partitioner=LabelInOut)]

build_experiments(experiment_params_list)

run_estimation_experiments(experiment_params_list)

graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=error, filename="LabelInOut")