
using Profile
include("../Experiments.jl")

# The goal of this file is to evaluate different estimation methods across different datasets.

datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
bounds_datasets = [human, aids, lubm80]

bounds_mix_scheme = [(Degree, 8), (QuasiStable, 8), (NeighborNodeLabels, 8), (NodeLabels, 8)]
mix_scheme = [(Degree, 8), (QuasiStable, 8), (NeighborNodeLabels, 8), (NodeLabels, 8)]

experiment_params = Vector{ExperimentParams}()
max_bounds_experiment_params = Vector{ExperimentParams}()
min_bounds_experiment_params = Vector{ExperimentParams}()
smaller_experiment_params = Vector{ExperimentParams}()
for dataset in bounds_datasets
    push!(max_bounds_experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=bounds_mix_scheme,
                                                max_cycle_size = -1,
                                                description = "COLOR (MaxMix32)"))

    push!(max_bounds_experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                max_cycle_size = -1,
                                                inference_max_paths = 10^30,
                                                summary_max_paths=1000,
                                                use_partial_sums = true,
                                                description = "BSK++"))
end
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                description = "COLOR \n(AvgMix32)"))
                                                
    push!(smaller_experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                description = "COLOR (AvgMix32)"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 1)],
                                                max_cycle_size = -1,
                                                description = "TradEst"))

    
end

println("Building...")

# build_experiments(experiment_params)
# build_experiments(max_bounds_experiment_params)

println("Estimating...")

# run_estimation_experiments(experiment_params; timeout=TIMEOUT_SEC)
# run_estimation_experiments(max_bounds_experiment_params; timeout=TIMEOUT_SEC)

comparison_methods =  ["alley", "alleyTPI", "wj", "impr", "jsub", "cs", "cset", "sumrdf", "lss"]
x_order = [string(data) for data in datasets]
bounds_x_order = [string(data) for data in bounds_datasets]
legend_order = [params.description for params in experiment_params][1:Int(length(experiment_params)/ length(datasets))]
max_bounds_legend_order = [params.description for params in max_bounds_experiment_params][1:Int(length(max_bounds_experiment_params)/ length(bounds_datasets))]
legend_order = vcat(legend_order, comparison_methods)

colors = [:red :yellow :maroon3 :palevioletred1 :dodgerblue :coral :palegreen :mediumpurple2 :darkgreen :cadetblue1 :goldenrod]

println("Graphing figures 3 and 4...")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                x_type = dataset,
                                                x_order = x_order,
                                                legend_order = legend_order,
                                                grouping=description,
                                                dimensions = (1550, 650),
                                                legend_pos=:outerright,
                                                legend_columns = 1,
                                                y_label="Relative Error log\$_{10}\$",
                                                group_colors = colors,
                                                filename="fig_3") # overall error

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
                                                y_type = runtime,
                                                x_type = dataset,
                                                x_order = x_order,
                                                legend_order = legend_order,
                                                grouping=description,
                                                dimensions = (1550, 650),
                                                legend_pos=:outerright,
                                                legend_columns = 1,
                                                y_label="Inference Latency log\$_{10}\$ (s)",
                                                group_colors = colors,
                                                filename="fig_4") # overall runtime

println("Graphing figures 5 and 6...")

graph_grouped_box_plot(max_bounds_experiment_params;
                                                ylims=[10^0, 10^30],
                                                y_ticks=[10^0, 10^5, 10^10, 10^15],
                                                y_type = estimate_error,
                                                x_type = dataset,
                                                x_order = bounds_x_order,
                                                legend_order = max_bounds_legend_order,
                                                grouping=description,
                                                dimensions = (600, 400),
                                                legend_pos=:topleft,
                                                legend_columns=1,
                                                # include_hline = false,
                                                y_label="Relative Error log\$_{10}\$",
                                                filename="fig_5") # bounds error

graph_grouped_box_plot(max_bounds_experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
                                                y_type = runtime,
                                                x_type = dataset,
                                                x_order = bounds_x_order,
                                                legend_order = max_bounds_legend_order,
                                                grouping=description,
                                                dimensions = (600, 400),
                                                legend_pos=:topright,
                                                y_label="Inference Latency log\$_{10}\$ (s)",
                                                filename="fig_6") # bounds runtime

comparison_methods =  ["alleyTPI", "sumrdf", "lss"]
x_order = [string(data) for data in datasets]
bar_legend_order = [params.description for params in smaller_experiment_params][1:Int(length(smaller_experiment_params)/ length(datasets))]
bar_legend_order = vcat(bar_legend_order, comparison_methods)
println("bar legend order: ", bar_legend_order)
bar_plot_colors = [:red :palevioletred1 :cadetblue1 :goldenrod]

println("Graphing figures 7 and 8")

graph_grouped_bar_plot(smaller_experiment_params;
                        grouping=description,
                        y_type=memory_footprint,
                        x_order = x_order,
                        legend_order = bar_legend_order,
                        ylims=[0, 10],
                        y_ticks = [1, 2, 3, 4, 5, 6, 7, 8],
                        legend_pos=:topleft,
                        dimensions = (900, 400),
                        scale_factor = 1000,
                        log_scale = true,
                        group_colors = bar_plot_colors,
                        y_label="Memory log\$_{10}\$ (KB)",
                        filename="fig_7") # overall memory

graph_grouped_bar_plot(smaller_experiment_params;
                        grouping=description,
                        y_type=build_time,
                        x_order = x_order,
                        legend_order = bar_legend_order,
                        legend_pos=:topleft,
                        ylims=[0, 11],
                        y_ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        dimensions = (900, 400),
                        scale_factor = 1000,
                        log_scale = true,
                        group_colors = bar_plot_colors,
                        y_label="Build Time log\$_{10}\$ (ms)",
                        filename="fig_8") # overall build time
