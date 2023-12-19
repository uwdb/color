using Plots.PlotMeasures
include("../Experiments.jl")

# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
data = [yeast]
max_cycles = 6
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, max_cycle_size=current_size) for current_dataset in data for current_size in 2:max_cycles]

# println("started building")
#build_experiments(experiment_params_list)
# println("started estimating")
#run_estimation_experiments(experiment_params_list)
println("started graphing")

graph_grouped_box_plot(experiment_params_list;
                        x_type=dataset,
                        ylims=[10^-20, 10^15],
                        y_ticks =[10^-20,10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        grouping=cycle_size,
                        y_type=estimate_error,
                        y_label="Relative Error 10^ (s)",
                        filename="cycles_size_error")

graph_grouped_box_plot(experiment_params_list;
                        x_type=dataset,
                        ylims=[10^-3, 10^1],
                        y_ticks =[10^-3, 10^-2, 10^-1, 1, 10],
                        grouping=cycle_size,
                        y_type=runtime,
                        y_label="Seconds 10^ (s)",
                        filename="cycles_size_runtime")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=memory_footprint,
                        ylims=[0, 15.5],
                        y_ticks =[3, 6, 9, 12, 15],
                        dimensions = (600, 400),
                        y_label="Memory Footprint (MB)",
                        grouping=cycle_size,
                        filename="cycles_size_memory")

graph_grouped_bar_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=build_time,
                        ylims = [0, 4],
                        y_ticks = [1, 2, 3, 4],
                        dimensions = (600, 400),
                        legend_pos = :topright,
                        y_label="Build Time (s)",
                        grouping=cycle_size,
                        filename="cycles_size_build_time")
