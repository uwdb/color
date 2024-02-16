using Plots.PlotMeasures
include("../Experiments.jl")

# The goal of this file is to evaluate the effect of using different values of maximum partial paths (varying the amount of sampling) during summary-building time, 
# across different datasets.

datasets::Vector{DATASET} = [yeast]
max_paths = [64, 256, 1024, 4096, 8192, 8192*4, 8192*4*4]

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, summary_max_paths=current_paths) for current_dataset in datasets for current_paths in max_paths]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=build_time,
                        ylims=[0, 200],
                        y_ticks=[50, 100, 150, 200],
                        y_label="Build Time",
                        grouping=summary_paths,
                        filename="summary-build-multi")

graph_grouped_box_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=estimate_error,
                        ylims=[10^-20, 10^15],
                        y_ticks =[10^-20, 10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        y_label="Estimate Error", grouping=summary_paths,
                        filename="summary-error-multi")

graph_grouped_box_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=runtime,
                        ylims=[.0001, 10],
                        y_ticks =[.0001, .001, .01, 1, 10],
                        y_label="Runtime",
                        grouping=summary_paths,
                        filename="summary-runtime-multi")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=memory_footprint,
                        ylims=[0, 30],
                        y_ticks =[5 ,10, 15, 20, 25, 30],
                        y_label="Memory Footprint",
                        grouping=summary_paths,
                        filename="summary-memory-multi")
