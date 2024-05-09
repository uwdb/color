using Plots.PlotMeasures
include("../Experiments.jl")

# The goal of this file is to evaluate summaries using different maximum stored cycle size statistics, on the same dataset.

datasets = [youtube]
max_cycles = 6
experiment_params_list = ExperimentParams[ExperimentParams(dataset=current_dataset, max_cycle_size=current_size) for current_dataset in datasets for current_size in 1:max_cycles]

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")

graph_grouped_box_plot(experiment_params_list;
                        x_type=cycle_size,
                        ylims=[10^-20, 10^15],
                        y_ticks =[10^-20,10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        grouping=dataset,
                        y_type=estimate_error,
                        dimensions = (600, 400),
                        legend_pos = nothing,
                        y_label="Relative Error log\$_{10}\$",
                        x_label="Maximum Cycle Size",
                        filename="fig_12") # cycles size error

graph_grouped_box_plot(experiment_params_list;
                        x_type=cycle_size,
                        ylims=[10^-3, 10^1],
                        y_ticks =[10^-3, 10^-2, 10^-1, 1, 10],
                        grouping=dataset,
                        y_type=runtime,
                        dimensions = (600, 400),
                        legend_pos = nothing,
                        y_label="Inference Latency log\$_{10}\$ (s)",
                        x_label="Maximum Cycle Size",
                        filename="cycles-size-runtime")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=cycle_size,
                        y_type=memory_footprint,
                        ylims=[0, 15.5],
                        y_ticks =[3, 6, 9, 12, 15],
                        dimensions = (600, 400),
                        y_label="Statistics Size (MB)",
                        x_label="Maximum Cycle Size",
                        grouping=dataset,
                        legend_pos = nothing,
                        filename="cycles-size-memory")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=cycle_size,
                        y_type=build_time,
                        ylims = [0, 4],
                        y_ticks = [1, 2, 3, 4],
                        dimensions = (600, 400),
                        y_label="Build Time (s)",
                        x_label="Maximum Cycle Size",
                        grouping=dataset,
                        legend_pos = nothing,
                        filename="cycles-size-build-time")
