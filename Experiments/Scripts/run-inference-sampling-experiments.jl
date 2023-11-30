using Plots.PlotMeasures
include("../Experiments.jl")

# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
datasets::Vector{DATASET} = [aids]
max_paths = 60
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, inference_max_paths=current_paths) for current_dataset in datasets for current_paths in 10:10:max_paths]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, x_label="dataset", y_label="accuracy", grouping=inference_paths, filename="inferencesampling")