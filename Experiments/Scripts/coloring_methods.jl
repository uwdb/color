
using Profile
include("../Experiments.jl")

datasets = [human, aids, yeast, hprd, dblp]
partitioners = [QuasiStable, Degree, Hash]
label_refining_rounds = [0, 1, 2, 3, 4]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    for partitioner in partitioners
        for refining_round in label_refining_rounds
            num_initial_partitions = Int(128/(2 ^ refining_round))
            push!(experiment_params, ExperimentParams(dataset=dataset, partitioner=partitioner,
                    num_colors = num_initial_partitions,
                    label_refining_rounds=refining_round, sampling_strategy=redistributive))
            println(num_initial_partitions)
        end
    end
end

build_experiments(experiment_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=technique, filename="coloring_strategies_comparison_equal_colors")

graph_grouped_box_plot(experiment_params; grouping=technique, y_type=runtime, filename="coloring_strategies_comparison_runtime_equal_colors")
