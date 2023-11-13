
using Profile
include("../Experiments.jl")

#datasets = [aids, yeast, hprd, dblp, youtube, wordnet]
datasets = [lubm80, aids, human]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(dataset=dataset))
    push!(experiment_params, ExperimentParams(dataset=dataset, num_colors=256))
end

#build_experiments(experiment_params)

run_estimation_experiments(experiment_params)

graph_grouped_boxplot_with_comparison_methods(experiment_params; ylims=[10^-5, 10^2],y_type = runtime, grouping=number_of_colors, y_label="Runtime (s)", filename="comparison_exps_runtime_2")
graph_grouped_boxplot_with_comparison_methods(experiment_params; y_type = estimate_error, grouping=number_of_colors, y_label="Relative Error", filename="comparison_exps_error_2")
