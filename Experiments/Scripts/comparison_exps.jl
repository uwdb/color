
using Profile
include("../Experiments.jl")

#datasets = [aids, yeast, hprd, dblp, youtube, wordnet]
datasets = [yeast]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=CorrDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "CorrQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=CorrDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 32), (NeighborNodeLabels, 32),(QuasiStable, 32), (NeighborNodeLabels, 32)],
                                                description = "CorrQ64N64"))
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
                                                description = "MaxH64"))
end

build_experiments(experiment_params)

run_estimation_experiments(experiment_params)

graph_grouped_boxplot_with_comparison_methods(experiment_params; ylims=[10^-5, 10^2], y_type = runtime, grouping=description, y_label="Runtime (s)", filename="comparison_exps_runtime_2")
graph_grouped_boxplot_with_comparison_methods(experiment_params; ylims=[10^-10, 10^15], y_type = estimate_error, grouping=description, y_label="Relative Error", filename="comparison_exps_error_2")
