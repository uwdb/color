using Plots.PlotMeasures
include("Experiments/Experiments.jl")

datasets::Vector{DATASET} = [aids]
# datasets::Vector{DATASET} = [aids, human, yeast, wordnet, youtube, dblp, patents]
# datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
max_cycles = 6

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_size) for current_dataset in datasets for current_size in 2:max_cycles]
println("started building")
for experiment_params in experiment_params_list
    build_times = [("Dataset", "Partitioner", "NumColors", "BuildTime", "MemoryFootprint")]
    dataset = experiment_params.dataset
    summary_params = experiment_params.summary_params
    data = load_dataset(dataset)
    summary_name = params_to_summary_filename(experiment_params)
    summary_file_location = "Experiments/SerializedSummaries/" * summary_name
    println("Building Color Summary: ", summary_name)
    # normal_results = @timed generate_color_summary(data, summary_params; verbose=1)
    # normal_results = @timed generate_color_summary(data, summary_params; verbose=0, detailed_cycles=false)
    results= @timed generate_color_summary(data, summary_params; verbose=1, detailed_cycles=false)
    # println("normal time: ", normal_results.time)
    println("detailed time: ", results.time)
    summary_size = Base.summarysize(results.value)
    serialize(summary_file_location, results.value)
    push!(build_times, (string(dataset),
                         string(summary_params.partitioner),
                         string(summary_params.num_colors),
                         string(results.time),
                         string(summary_size)))
    results_filename = params_to_results_filename(experiment_params)
    result_file_location = "Experiments/Results/Build_" * results_filename
    writedlm(result_file_location, build_times, ",")
end
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=cycle_size, filename="detailed-sample-experiment")