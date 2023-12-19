using Plots.PlotMeasures
include("../Experiments.jl")

# datasets::Vector{DATASET} = [aids, human, lubm80, yago, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
data::Vector{DATASET} = [yeast]
max_paths = [10, 100, 1000, 10000]
experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, inference_max_paths=current_paths) for current_dataset in data for current_paths in max_paths]
println("started building")
#build_experiments(experiment_params_list)
println("started estimating")
#run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=estimate_error,
                        ylims=[10^-20, 10^20],
                        y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        dimensions = (600, 550),
                        legend_pos = :topleft,
                        x_label="Maximum Inference Paths",
                        y_label="Estimate Error 10^",
                        grouping=inference_paths,
                        filename="inference-paths-error")

graph_grouped_box_plot(experiment_params_list,
                        x_type=dataset,
                        y_type=runtime,
                        ylims=[.0001, 100],
                        y_ticks=[.001, .01, .1, 1, 10, 100],
                        dimensions = (600, 550),
                        legend_pos = :topleft,
                        x_label="Maximum Inference Paths",
                        y_label="Runtime 10^ (s)",
                        grouping=inference_paths,
                        filename="inference-paths-runtime")
