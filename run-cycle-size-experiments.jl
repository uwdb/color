# TODO:
# - turn sampling up/down (build and inference)
# - cycle length effects
# - reverify ground truth wrt G-captures
# - rerun initial g-care benchmarks and verify old results
# - recreate initial G-Care Benchmark results

# My tasks:
# - turn sampling up/down
# - cycle length effects


include("Experiments/build_color_summaries.jl")
include("Experiments/get_true_cardinalities.jl")
include("Experiments/load_datasets.jl")
include("Experiments/load_querysets.jl")
include("Experiments/run_estimators.jl")
include("Experiments/graph_results.jl")
include("Experiments/utils.jl")

datasets::Vector{DATASET} = [aids, wordnet, lubm80, human]
max_cycle_size = 5
# experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_size) for current_dataset in datasets for current_size in 0:max_cycle_size]

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=aids, partitioner=QuasiStable, max_cycle_size=6),
                                                    ExperimentParams(dataset=aids, partitioner=QuasiStable, max_cycle_size=5),
                                                    ExperimentParams(dataset=aids, partitioner=QuasiStable, max_cycle_size=4),
                                                    ExperimentParams(dataset=aids, partitioner=QuasiStable, max_cycle_size=3),
                                                    ExperimentParams(dataset=aids, partitioner=QuasiStable, max_cycle_size=2)]

                                                    # ExperimentParams(dataset=wordnet, partitioner=QuasiStable, max_cycle_size=6),
                                                    # ExperimentParams(dataset=wordnet, partitioner=QuasiStable, max_cycle_size=5),
                                                    # ExperimentParams(dataset=wordnet, partitioner=QuasiStable, max_cycle_size=4),
                                                    # ExperimentParams(dataset=wordnet, partitioner=QuasiStable, max_cycle_size=3),
                                                    # ExperimentParams(dataset=wordnet, partitioner=QuasiStable, max_cycle_size=2),

                                                    # ExperimentParams(dataset=lubm80, partitioner=QuasiStable, max_cycle_size=6),
                                                    # ExperimentParams(dataset=lubm80, partitioner=QuasiStable, max_cycle_size=5),
                                                    # ExperimentParams(dataset=lubm80, partitioner=QuasiStable, max_cycle_size=4),
                                                    # ExperimentParams(dataset=lubm80, partitioner=QuasiStable, max_cycle_size=3),
                                                    # ExperimentParams(dataset=lubm80, partitioner=QuasiStable, max_cycle_size=2),

                                                    # ExperimentParams(dataset=human, partitioner=QuasiStable, max_cycle_size=6),
                                                    # ExperimentParams(dataset=human, partitioner=QuasiStable, max_cycle_size=5),
                                                    # ExperimentParams(dataset=human, partitioner=QuasiStable, max_cycle_size=4),
                                                    # ExperimentParams(dataset=human, partitioner=QuasiStable, max_cycle_size=3),
                                                    # ExperimentParams(dataset=human, partitioner=QuasiStable, max_cycle_size=2)]
build_experiments(experiment_params_list)

run_estimation_experiments(experiment_params_list)

graph_grouped_box_plots_cycles(experiment_params_list, graph_title="cyclesize")

# tmux new -A -s experiments