using Graphs: DiGraph
using StatsPlots
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!
using Test
using Graphs
using AutoHashEquals
include("../src/datasets.jl")
include("../src/QuasiStableCardinalityEstimator.jl")

first_bool_path::BoolPath = [true, true]
copy_bool_path::BoolPath = [true, true]
first_colors::StartEndColorPair = [1, 1]
copy_colors::StartEndColorPair = [1, 1]

first_cycle_info = CyclePathAndColors(first_bool_path, first_colors)
copy_cycle_info = CyclePathAndColors(copy_bool_path, copy_colors)

test_dict::Dict{CyclePathAndColors, Float64} = Dict()
test_dict[first_cycle_info] = 1
test_dict[copy_cycle_info] = 1
println(test_dict)
println(isequal(first_cycle_info, copy_cycle_info))


first_string = "hi"
copy_string = "hi"

struct StringStorer
    a::String
end

first_struct = StringStorer(first_string)
second_struct = StringStorer(copy_string)
string_dict::Dict{StringStorer, Float64} = Dict()
string_dict[first_struct] = 1
string_dict[second_struct] = 1
println(string_dict)
