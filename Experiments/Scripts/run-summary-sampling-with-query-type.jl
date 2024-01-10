using Plots.PlotMeasures
include("../Experiments.jl")

current_dataset = yeast
max_paths = 300

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, summary_max_paths=current_paths) for current_paths in 10:30:max_paths]
println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=query_type, y_type=estimate_error, grouping=summary_paths, filename="summarysamplesquerytypesyeast")