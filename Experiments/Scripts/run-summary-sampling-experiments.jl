using Plots.PlotMeasures
include("../Experiments.jl")


datasets::Vector{DATASET} = [aids, wordnet, lubm80, human]
max_paths = 1000
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, summary_max_paths=current_paths) for current_dataset in datasets for current_paths in 0:200:max_paths]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=error, grouping=summary_paths, filename="summarysamples")
