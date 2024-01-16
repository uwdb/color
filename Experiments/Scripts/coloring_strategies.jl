
using Profile
include("../Experiments.jl")

datasets = [human, eu2005, dblp, youtube]
datasets = [human, hprd]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 16), (NeighborNodeLabels, 8), (NodeLabels, 8)],
                                                description = "AvgMix32"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 32), (NeighborNodeLabels, 16), (NodeLabels, 16)],
                                                description = "AvgMix64"))
#=
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64), (NeighborNodeLabels, 32), (NodeLabels, 32)],
                                                description = "AvgMix128"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(NodeLabels, 64)],
                                                description = "AvgN64"))

   push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(NeighborNodeLabels, 64)],
                                                description = "AvgNNL64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Degree, 64)],
                                                description = "AvgD64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "AvgQ64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                description = "AvgH64")) =#
end

build_experiments(experiment_params)

run_estimation_experiments(experiment_params; timeout=TIMEOUT_SEC)

x_order = [string(data) for data in datasets]
legend_order = [params.description for params in experiment_params][1:Int(length(experiment_params)/length(datasets))]

graph_grouped_box_plot(experiment_params;
                        ylims=[10^-5, 10^4],
                        y_ticks=[10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
                        y_type = runtime,
                        x_type = dataset,
                        x_order = x_order,
                        legend_order=legend_order,
                        grouping=description,
                        dimensions = (700, 450),
                        legend_pos=:top,
                        y_label="Inference Latency 10^ (s)",
                        filename="colorings_runtime")

graph_grouped_box_plot(experiment_params;
                        ylims=[10^-21, 10^21],
                        y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                        y_type = estimate_error,
                        x_type = dataset,
                        x_order = x_order,
                        legend_order=legend_order,
                        grouping=description,
                        dimensions = (700, 450),
                        legend_pos=:topleft,
                        y_label="Relative Error 10^",
                        filename="colorings_error")


graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=memory_footprint,
                        x_order = x_order,
                        legend_order=legend_order,
                        ylims=[0, 50],
                        y_ticks = [10, 20, 30, 40, 50],
                        legend_pos=:topright,
                        dimensions = (1000, 550),
                        y_label="Memory (MBs)",
                        filename="colorings_memory")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=build_time,
                        x_order = x_order,
                        legend_order=legend_order,
                        ylims=[0, 1600],
                        y_ticks = [200, 400, 600, 800, 1000, 1200, 1400, 1600],
                        dimensions = (1000, 550),
                        y_label="Build Time (s)",
                        filename="colorings_build_time")
