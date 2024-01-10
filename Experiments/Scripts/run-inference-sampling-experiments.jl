using Plots.PlotMeasures
include("../Experiments.jl")

datasets = [yeast]
max_paths = [10, 100, 500, 2000, 10000]
experiment_params_list = ExperimentParams[]
for dataset in datasets
    for current_paths in max_paths
        push!(experiment_params_list, ExperimentParams(dataset=dataset, inference_max_paths=current_paths, sampling_strategy = redistributive, description="Importance"))
        push!(experiment_params_list, ExperimentParams(dataset=dataset, inference_max_paths=current_paths, sampling_strategy = uniform, description="Uniform"))
    end
end

println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list,
                        x_type=inference_paths,
                        y_type=estimate_error,
                        ylims=[10^-20, 10^20],
                        y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        dimensions = (600, 400),
                        legend_pos = :topleft,
                        x_label="Maximum Inference Paths",
                        y_label="Estimate Error 10^",
                        grouping=description,
                        filename="inference-paths-error")

graph_grouped_box_plot(experiment_params_list,
                        x_type=inference_paths,
                        y_type=runtime,
                        ylims=[.0001, 100],
                        y_ticks=[.001, .01, .1, 1, 10, 100],
                        dimensions = (600, 400),
                        legend_pos = :topleft,
                        x_label="Maximum Inference Paths",
                        y_label="Runtime 10^ (s)",
                        grouping=description,
                        filename="inference-paths-runtime")
