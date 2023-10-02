@enum DATASET  aids human lubm80 yago yeast hprd wordnet dblp youtube eu2005 patents

const IS_GCARE_DATASET = Dict(aids=>true, human=>true, lubm80=>true, yago=>true,
     yeast=>false, hprd=>false, wordnet=>false, dblp=>false, youtube=>false, eu2005=>false, patents=>false)

struct ExperimentParams
    dataset::DATASET
    summary_params::ColorSummaryParams
    inference_max_paths::Int
    use_partial_sums::Bool

    function ExperimentParams(;dataset::DATASET,  num_colors::Int=64, max_cycle_size=4, summary_max_paths=1000,
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
    name *= params_to_string(experiment_params.summary_params) *  ".obj"
    return name
end
