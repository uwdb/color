
using Profile
include("../Experiments.jl")

datasets = [human, aids, yeast, hprd, dblp]
partitioners = [QuasiStable]
label_refining_rounds = [0, 2, 4, 6]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    for partitioner in partitioners
        for refining_round in label_refining_rounds
            push!(experiment_params, ExperimentParams(dataset=dataset, partitioner=partitioner,
                    label_refining_rounds=refining_round, inference_max_paths=1024))
        end
    end
end

#build_experiments(experiment_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=technique, filename="coloring_strategies_comparison_best_2")
