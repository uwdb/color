
using Profile
include("../Experiments.jl")

#datasets = [human, aids, yeast, hprd, dblp, wordnet]
datasets = [aids]

experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
num_colors = 128
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset,  num_colors=num_colors))
#    for sampling_strategy in instances(SAMPLING_STRATEGY)
    for sampling_strategy in [redistributive, loop_vec, online]
            push!(experiment_params, ExperimentParams(dataset=dataset, num_colors=num_colors, sampling_strategy=sampling_strategy, inference_max_paths=512))
    end
end

#build_experiments(build_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=sampling_type, filename="sampling_strategies_comparison_error")

graph_grouped_box_plot(experiment_params; y_type = runtime, grouping=sampling_type, filename="sampling_strategies_comparison_runtime")
