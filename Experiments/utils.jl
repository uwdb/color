@enum DATASET  aids human lubm80 yago yeast hprd wordnet dblp youtube eu2005 patents

const IS_GCARE_DATASET = Dict(aids=>true, human=>true, lubm80=>true, yago=>true,
     yeast=>false, hprd=>false, wordnet=>false, dblp=>false, youtube=>false, eu2005=>false, patents=>false)

struct ExperimentParams
    dataset::DATASET
    summary_params::ColorSummaryParams
    inference_max_paths::Int
    # Whether to consider just the shortest path when calculating cycle probabilities or
    # all simple paths.
    only_shortest_path_cycle::Bool
    use_partial_sums::Bool
    sampling_strategy::SAMPLING_STRATEGY
    description::String
    n_replications::Int

    function ExperimentParams(;dataset::DATASET, max_cycle_size=6,
        only_shortest_path_cycle=false, summary_max_paths=50000,
        partitioning_scheme::Vector{Tuple{PARTITIONER, Int}} = [(QuasiStable, 64)], weighting=true,
        inference_max_paths=500, use_partial_sums=true,
        sampling_strategy=redistributive, proportion_updated=0.0, proportion_deleted=0.0,
        deg_stats_type::Type=AvgDegStats, description="", n_replications=3)
        return new(dataset,
                    ColorSummaryParams(deg_stats_type=deg_stats_type,
                                        max_cycle_size=max_cycle_size,
                                        max_partial_paths=summary_max_paths,
                                        partitioning_scheme=partitioning_scheme,
                                        weighting=weighting,
                                        proportion_updated=proportion_updated,
                                                       proportion_deleted=proportion_deleted),
                    inference_max_paths,
                    only_shortest_path_cycle,
                    use_partial_sums,
                    sampling_strategy,
                    description,
                    n_replications
               )
    end
end

function params_to_results_filename(experiment_params::ExperimentParams)
    name = string(experiment_params.dataset) * "_"
    name *= params_to_string(experiment_params.summary_params) * "_"
    name *= string(experiment_params.inference_max_paths) * "_"
    name *= string(experiment_params.only_shortest_path_cycle) * "_"
    name *= string(experiment_params.use_partial_sums) * "_"
    # name *= string(experiment_params)
    name *= string(experiment_params.sampling_strategy) * ".csv"
    return name
end

function params_to_summary_filename(experiment_params::ExperimentParams)
    name = string(experiment_params.dataset) * "_"
    name *= params_to_string(experiment_params.summary_params) *  ".obj"
    return name
end
