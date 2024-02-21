using qColoringCardinality
using Distributions
using DataStructures: counter, Dict, Set, Vector, inc!
using Test
using Graphs
include("../src/CardinalityWithColors.jl")


# run subtests
include("graphsize.jl")
include("statcheck.jl")
