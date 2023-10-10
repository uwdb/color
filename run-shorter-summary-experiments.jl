using Plots.PlotMeasures 
include("Experiments/build_color_summaries.jl")
include("Experiments/get_true_cardinalities.jl")
include("Experiments/load_datasets.jl")
include("Experiments/load_querysets.jl")
include("Experiments/run_estimators.jl")
include("Experiments/graph_results.jl")
include("Experiments/utils.jl")

# datasets::Vector{DATASET} = [aids, wordnet, lubm80, human]
# max_partial_paths = 10000

# experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=2),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=10),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=50),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=250),
#                                                     ExperimentParams(dataset=aids, partitioner=QuasiStable, inference_max_paths=1250)]

datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, youtube]
max_paths = 60

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, summary_max_paths=current_paths) for current_dataset in datasets for current_paths in 2:10:max_paths]

build_experiments(experiment_params_list)

run_estimation_experiments(experiment_params_list)

graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=error, grouping=summary_paths, filename="shortersummarysamples")