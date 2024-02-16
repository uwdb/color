#!/bin/bash

# This script runs the experiments included in the submitted paper for COLOR and saves the associated figures.

# Generate figures 2, 3, 4, 5, 6
julia Experiments/Scripts/comparison_exps.jl

# Generate figure 7
julia Experiments/Scripts/coloring_strategies.jl

# Generate figure 8
julia Experiments/Scripts/construction_scaling.jl

# Generate figure 9
julia Experiments/Scripts/proportion_updated.jl

# Generate figure 10
julia Experiments/Scripts/max_inference_paths.jl

# Generate figure 11
julia Experiments/Scripts/max_cycle_size.jl

# Generate figure 12
julia Experiments/Scripts/query_path_width_build.jl