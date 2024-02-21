
using Profile
include("../Experiments.jl")

# The goal of this file is to compare the runtime and accuracy for summaries using different partioning schemes on the same dataset.

datasets = [aids]
partitioning_schemes = [
                        [(Degree, 64)],
                        [(NeighborNodeLabels, 64)],
                        [(QuasiStable, 64)],
                        [(QuasiStable, 32), (NeighborNodeLabels, 32)],
                        [(Hash, 64)],
                        [(Degree, 8), (QuasiStable, 32), (NeighborNodeLabels, 24)],
                        [(Degree, 8), (NeighborNodeLabels, 24), (QuasiStable, 32)]]
experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    for scheme in partitioning_schemes
            push!(experiment_params, ExperimentParams(dataset=dataset, partitioning_scheme=scheme))
    end
end

build_experiments(experiment_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=technique, filename="coloring_strategies_comparison_equal_colors")

graph_grouped_box_plot(experiment_params; grouping=technique, y_type=runtime, filename="coloring_strategies_comparison_runtime_equal_colors")
