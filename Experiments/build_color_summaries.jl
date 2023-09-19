using Serialization: serialize
using DelimitedFiles: writedlm
include("load_datasets.jl")
include("utils.jl")

function build_experiments(experiment_params_list::Vector{ExperimentParams})
    build_times = [("Dataset", "Partitioner", "NumColors", "BuildTime")]
    for experiment_params in experiment_params_list
        dataset = experiment_params.dataset
        summary_params = experiment_params.summary_params
        data = load_dataset(dataset)
        summary_name = string(dataset) * "_" * params_to_string(params)
        summary_file_location = "Experiments/SerializedSummaries/" * summary_name * ".obj"
        println("Building Color Summary: ", summary_name)
        results = @timed generate_color_summary(data, summary_params; verbose=1)
        serialize(summary_file_location, results.value)
        push!(build_times, (string(dataset), string(summary_params.partitioner),
                             string(summary_params.num_colors), string(results.time)))
        writedlm("Experiments/Results/build_times.csv", build_times, ",") #TODO: use the params file name
    end
end
