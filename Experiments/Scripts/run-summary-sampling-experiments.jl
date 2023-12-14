using Plots.PlotMeasures
include("../Experiments.jl")


# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
datasets::Vector{DATASET} = [aids, human, dblp]
# max_paths = [32, 64, 128, 256, 512, 1024, 2048, 20000, 100000]
max_paths = [32, 64, 128, 256, 512, 1024, 2048]

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, summary_max_paths=current_paths) for current_dataset in datasets for current_paths in max_paths]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
# graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, x_label="dataset", y_label="accuracy", grouping=summary_paths, filename="summarysamples")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=build_time, y_lims=[0, 100], x_label="Maximum Paths During Summary Building", y_label="Build Time", grouping=summary_paths, filename="summary-build-multi")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, x_label="Maximum Paths During Summary Building", y_label="Estimate Error", grouping=summary_paths, filename="summary-error-multi")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=runtime, y_lims=[0, 0.5], x_label="Maximum Paths During Summary Building", y_label="Runtime", grouping=summary_paths, filename="summary-runtime-multi")
graph_grouped_bar_plot(experiment_params_list, x_type=dataset, y_type=memory_footprint, y_lims=[0, 20], x_label="Maximum Paths During Summary Building", y_label="Memory Footprint", grouping=summary_paths, filename="summary-memory-multi")