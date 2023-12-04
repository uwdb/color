using Plots.PlotMeasures
include("../Experiments.jl")

# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
datasets::Vector{DATASET} = [aids]
max_cycles = 6
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_size) for current_dataset in datasets for current_size in 2:max_cycles]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_bar_plot(experiment_params_list, x_type=cycle_size, y_type=build_time, x_label="Size of Cycles for Stored Probability", y_label="Build Time", grouping=dataset, filename="cycle-build-aids")
graph_grouped_box_plot(experiment_params_list, x_type=cycle_size, y_type=estimate_error, x_label="Size of Cycles for Stored Probability", y_label="Estimate Error", grouping=dataset, filename="cycle-error-aids")
graph_grouped_bar_plot(experiment_params_list, x_type=cycle_size, y_type=runtime, y_lims=[0, 0.02], x_label="Size of Cycles for Stored Probability", y_label="Runtime", grouping=dataset, filename="cycle-runtime-aids")
graph_grouped_bar_plot(experiment_params_list, x_type=cycle_size, y_type=memory_footprint, y_lims=[0, 4], x_label="Size of Cycles for Stored Probability", y_label="Memory Footprint", grouping=dataset, filename="cycle-memory-aids")