# TODO:
# - turn sampling up/down (build and inference)
# - cycle length effects
# - reverify ground truth wrt G-captures
# - rerun initial g-care benchmarks and verify old results
# - recreate initial G-Care Benchmark results

# My tasks:
# - turn sampling up/down
# - cycle length effects

using Plots.PlotMeasures
include("Experiments/Experiments.jl")

# datasets::Vector{DATASET} = [aids, wordnet, lubm80, human]
# max_partial_paths = 10000

# experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=2),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=10),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=50),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=250),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=1250)]

datasets::Vector{DATASET} = [aids, human, yeast, wordnet, youtube, dblp, patents]
max_paths = 60

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, inference_max_paths=current_paths) for current_dataset in datasets for current_paths in 2:10:max_paths]
println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=inference_paths, filename="inferencesampling")