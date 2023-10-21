using Plots.PlotMeasures
include("Experiments/Experiments.jl")

# datasets::Vector{DATASET} = [aids, wordnet, lubm80, human]
# max_partial_paths = 10000

# experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=2),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=10),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=50),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=250),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=1250)]

# datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
datasets::Vector{DATASET} = [aids]

max_paths = 1000

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, summary_max_paths=current_paths) for current_dataset in datasets for current_paths in 250:250:max_paths]
println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=summary_paths, filename="summarysamples")