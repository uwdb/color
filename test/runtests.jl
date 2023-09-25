using qColoringCardinality
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!
using Test
using Graphs
include("../Source/CardinalityWithColors.jl")


# run subtests
include("graphsize.jl")
include("statcheck.jl")
