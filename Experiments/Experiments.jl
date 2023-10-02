# Overall Experiments Harness Include File
using Serialization: serialize, deserialize
using StatsPlots
using CSV, DataFrames
using DelimitedFiles: writedlm

include(pwd() * "/Source/CardinalityWithColors.jl")
include(pwd() * "/Experiments/utils.jl")
include(pwd() * "/Experiments/load_datasets.jl")
include(pwd() * "/Experiments/load_querysets.jl")
include(pwd() * "/Experiments/build_color_summaries.jl")
include(pwd() * "/Experiments/get_true_cardinalities.jl")
include(pwd() * "/Experiments/run_estimators.jl")
include(pwd() * "/Experiments/graph_results.jl")
