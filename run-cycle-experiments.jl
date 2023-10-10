using Plots.PlotMeasures 
include("Experiments/build_color_summaries.jl")
include("Experiments/get_true_cardinalities.jl")
include("Experiments/load_datasets.jl")
include("Experiments/load_querysets.jl")
include("Experiments/run_estimators.jl")
include("Experiments/graph_results.jl")
include("Experiments/utils.jl")

# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
max_cycles = 6

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_size) for current_dataset in datasets for current_size in 2:max_cycles]

build_experiments(experiment_params_list)

run_estimation_experiments(experiment_params_list)

graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=error, grouping=cycle_size, filename="justaidsagain")