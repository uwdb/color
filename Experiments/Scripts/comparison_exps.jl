
using Profile
include("../Experiments.jl")

#datasets = [human, aids, lubm80, yeast, hprd, dblp, youtube, eu2005, patents, wordnet]
datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
#datasets = [aids]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "AvgQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 32), (NeighborNodeLabels, 32),(QuasiStable, 32), (NeighborNodeLabels, 32)],
                                                description = "AvgQ64N64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=MinDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "MinQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "MaxQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                description = "BSK"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 1)],
                                                max_cycle_size = -1,
                                                description = "IndEst"))
end

build_experiments(experiment_params)

#run_estimation_experiments(experiment_params)

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-2, 10^0, 10^2, 10^4],
                                                y_type = runtime,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:top,
                                                y_label="Runtime (10^ s)",
                                                filename="runtime")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:bottomleft,
                                                y_label="Relative Error (10^)",
                                                filename="error")


graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                x_type = query_size,
                                                y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:bottomleft,
                                                y_label="Relative Error (10^)",
                                                filename="error-query-size")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=memory_footprint,
                        ylims=[0, 30],
                        y_ticks = [5, 10, 15, 20, 25, 30],
                        dimensions = (1000, 550),
                        y_label="Memory (MBs)",
                        filename="memory")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=build_time,
                        ylims=[0, 100],
                        y_ticks = [20, 40, 60, 80, 100],
                        dimensions = (1000, 550),
                        y_label="Build Time (s)",
                        filename="build_time")
