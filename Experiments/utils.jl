include("load_datasets.jl")


struct ExperimentParams
    dataset::DATASET
    summary_params::ColorSummaryParams
    inference_max_paths::Int
    use_partial_sums::Bool

    function ExperimentParams(;dataset::Dataset,  num_colors::Int=64, max_cycle_size=4, summary_max_paths=1000,
        partitioner::PARTITIONER = QuasiStable, weighting=true, inference_max_paths=500, use_partial_sums=true)
        return new(dataset, ColorSummaryParams(num_colors=num_colors,
                                                       max_cycle_size=max_cycle_size,
                                                       max_partial_paths=summary_max_paths,
                                                       partitioner=partitioner,
                                                       weighting=weighting),
                    inference_max_paths,
                    use_partial_sums
               )
    end
end

function params_to_results_filename(experiment_params::ExperimentParams)
    name = string(experiment_params.dataset) * "_"
    name *= string(experiment_params.summary_params) * "_"
    name *= string(experiment_params.inference_max_paths) * "_"
    name *= string(experiment_params.use_partial_sums) * ".csv"
    return name
end

function params_to_summary_filename(experiment_params::ExperimentParams)
    name = string(experiment_params.dataset) * "_"
    name *= string(experiment_params.summary_params) *  ".obj"
    return name
end
