
using Profile
include("../Experiments.jl")

# The goal of this file is to evaluate the effect of using sampling strategies, across different datasets.

datasets = [human, aids, yeast, hprd, dblp, wordnet]
datasets = [youtube]
experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset))
    for sampling_strategy in instances(SAMPLING_STRATEGY)
        push!(experiment_params, ExperimentParams(dataset=dataset, sampling_strategy=sampling_strategy, inference_max_paths=500))
    end
end

#build_experiments(build_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=sampling_type, filename="sampling_strategies_comparison")
