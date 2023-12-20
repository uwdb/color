# Overall Experiments Harness Include File
using Serialization: serialize, deserialize
using Plots
using Plots.PlotMeasures
using StatsPlots
using CSV, DataFrames
using Parquet2: Dataset
using DelimitedFiles: writedlm
using BenchmarkTools
using Random
using Printf
using SharedArrays
using WeakRefStrings
using Distributed

include("../Source/CardinalityWithColors.jl")
include("utils.jl")
include("load_datasets.jl")
include("load_querysets.jl")
include("build_color_summaries.jl")
include("get_true_cardinalities.jl")
include("run_estimators.jl")
include("graph_results.jl")
@everywhere include("../Source/CardinalityWithColors.jl")
@everywhere include("utils.jl")
@everywhere include("load_datasets.jl")
@everywhere include("build_color_summaries.jl")
@everywhere include("run_estimators.jl")
@everywhere using SharedArrays
@everywhere using WeakRefStrings
@everywhere using DelimitedFiles: writedlm
@everywhere using Parquet2: Dataset
@everywhere using Random
@everywhere using CSV, DataFrames
@everywhere using Serialization: serialize, deserialize
