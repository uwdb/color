# COLOR
<!-- # qColoringCardinality.jl -->

<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
<!-- ![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
<!-- [![build](https://github.com/mkyl/qColoringCardinality.jl/workflows/CI/badge.svg)](https://github.com/mkyl/qColoringCardinality.jl/actions?query=workflow%3ACI) -->
<!-- travis-ci.com badge, uncomment or delete as needed, depending on whether you are using that service. -->
<!-- [![Build Status](https://travis-ci.com/mkyl/qColoringCardinality.jl.svg?branch=master)](https://travis-ci.com/mkyl/qColoringCardinality.jl) -->
<!-- Coverage badge on codecov.io, which is used by default. -->
[![codecov.io](http://codecov.io/github/mkyl/qColoringCardinality.jl/coverage.svg?branch=master)](http://codecov.io/github/mkyl/qColoringCardinality.jl?branch=master)
<!-- Documentation -- uncomment or delete as needed -->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mkyl.github.io/qColoringCardinality.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://mkyl.github.io/qColoringCardinality.jl/dev)
-->

COLOR is a framework for improving traditional cardinality estimators in graph databases by applying graph colorings to produce a compact summary of the data graph. The project uses the Julia language for its implementation and provides tools to build graphs and run estimations.

## Installation

To set up a local copy of the project, complete the following steps:

1. Download the repository.
```
$ git clone
```

2. Include the necessary Julia Packages.
```
$ julia
julia> using Pkg;
julia> Pkg.instantiate();
```

3. Download query graphs, data graphs, and true cardinalities from [G-Care](https://github.com/yspark-dblab/gcare) and [In-Memory Subgraph Matching](https://github.com/RapidsAtHKUST/SubgraphMatching), also available as [zipped files](https://drive.google.com/drive/folders/1pjJz9ahXFEd3Nd1OxqLA2YNnXGuCVpEp?usp=sharing).

## API

Generally, the steps of a cardinality estimation will be:
1. Create a `QueryGraph`
2. Create a `DataGraph`
3. Create a `ColorSummary` of the `DataGraph`
4. Perform a cardinality estimation using the `QueryGraph` and the `ColorSummary`
5. Find the exact cardinality using the `DataGraph` and `QueryGraph` to evaluate the estimation.

### Graph Creation
We define `PropertyGraphs` which contain information about edge and vertex labels. `PropertyGraphs` have two types: `QueryGraphs` and `DataGraphs`. 

To create a `QueryGraph`, an initial graph with default wildcard (-1) labels can be constructed using either the number of vertices or an existing `DiGraph`, then the labels can be updated after construction.
```julia
# initialize a QueryGraph using the number of vertices...
q = QueryGraph(3)
# ... or initialize with an existing DiGraph
q = QueryGraph(DiGraph(3))
# now update the labels using functions from PropertyGraph.jl
add_labeled_edge!(q, (1, 2), 2)
update_data_labels!(q, 1, 3)
update_node_labels!(q, 2, [3,4])
```

`DataGraph`s are also created the same way with default wildcard (-1) labels. The main difference for a `DataGraph` is that its vertex data labels are not arbitrary and cannot be changed - instead, data labels are equivalent to `vertex_id - 1`.
```julia
# initialize a DataGraph using the number of vertices...
d = DataGraph(3)
# ... or initialize with an existing DiGraph
d = DataGraph(DiGraph(3))
# now update the labels using functions from PropertyGraph.jl
add_labeled_edge!(d, (1, 2), 2)
update_node_labels!(d, 2, [3,4])
```

Functions to convert graph files formatted like [G-Care](https://github.com/yspark-dblab/gcare) or [In-Memory Subgraph Matching](https://github.com/RapidsAtHKUST/SubgraphMatching) are also included in `src/Datasets.jl`: 
```julia
# convert a .txt file (following G-Care format) into a PropertyGraph
d1 = load_dataset("data.txt")
q1 = load_query("query.txt")
# convert a .graph file (following Subgraph-Matching format) into a PropertyGraph
d2 = load_dataset("data.graph",subgraph_matching_data=true)
q2 = load_query("query.graph", subgraph_matching_data=true)
```

### Summary Building

A lifted `ColorSummary` describing the overall `DataGraph` is necessary before performing cardinality estimation for a specified `QueryGraph`. To do this, simply use the `generate_color_summary` function from the `ColorSummary.jl` file:
```julia
d = load_dataset("data.txt")
summary = generate_color_summary(d)
```
Optional parameters for the summary-building can be changed such as the type of coloring to use or the amount of sampling to allow, but more information can be found in the code documentation.

### Cardinality Estimation

After obtaining a lifted `ColorSummary`, a cardinality estimation can be achieved for any `QueryGraph` by using the `get_cardinality_bounds` function from the `QuasiStableCardinalityEstimator.jl` file:
```julia
d = load_dataset("data.txt")
q = load_queryset("query.txt")
summary = generate_color_summary(d)
estimate = get_cardinality_bounds(q, summary)
```
The resulting estimate will be a singular `Float64` value.

Optional parameters for the estimation can be changed such as the sampling strategy or how to handle cycles, but more information can be found in the code documentation.


### Exact Cardinality
To calculate the accuracy of an estimation, we compare the result to the exact cardinality of the QueryGraph in the DataGraph. To do this, a `get_exact_size` function is included in the `ExactSizeCalculator.jl` file:
```julia
d = load_dataset("data.txt")
q = load_queryset("query.txt")
exact_cardinality = get_exact_size(q, d)
```
The resulting cardinality will be a singular `Float64` value.

Optional parameters for the calculation can be changed such as the timeout, but more information can be found in the code documentation.

The `load_true_cardinality` function from `src/Datasets.jl` also obtains the cardinality results from a given cardinality file (where the only element in the file is a singular cardinality value):
```julia
exact_cardinality = load_true_cardinality("cardinality.txt")
```

## Scripts

There are a variety of Julia scripts provided which perform different experiments then save figures presenting the results to the `Experiments/Results/Figures` folder. These scripts are stored in the `Experiments/Scripts` folder.

For example, to find how using different maximum cycle sizes in the summary affects the cardinality estimation, the `max_cycle_size.jl` script can be called from the main directory:
```
$ julia Experiments/Scripts/max_cycle_size.jl
```

The bash script `run_submitted_experiments.sh` in the `Experiments` folder is also included and will run the experiments described in the submitted paper. The script will run all included experiments then save all the corresponding figures. Any figures that are included in the paper have their file name match their figure number (i.e. the figure presenting the effect of inference sampling on relative error is named `fig_10.png`). The experiments included in this bash script are:
- `degree_variance_exps.jl` (figure 2)
- `comparison_exps.jl` (figures 3, 4, 5, 6, 7, 8)
- `coloring_strategies.jl` (figure 9)
- `construction_scaling.jl` (figure 10)
- `proportion_updated.jl` (figure 11)
- `max_cycle_size.jl` (figure 12)
- `query_path_width_build.jl` (figure 13)
- `max_inference_paths.jl` (figure 14)

The script can be called from the main directory:
```
$ Experiments/run_submitted_experiments.sh
```
