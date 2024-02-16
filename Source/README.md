# Source Files

The following files included in this folder implement different parts of the framework:

- `CardinalityWithColors.jl` - This file contains definitions for types that are used throughout the implementation, defines `ColorSummaryParams` used when building a `ColorSummary`, and includes all other files necessary for the cardinality estimator.

- `ColoringMethods.jl` - This file contains implementations of different ways to color the nodes in a graph.

- `ColorSummary.jl` - This file defines the `ColorSummary` struct used for the cardinality estimation and includes functions necessary to create one from a `DataGraph`.

- `Datasets.jl` - This file defines functions that can be used to obtain a `DataGraph`, `Querygraph`, or exact cardinality from a given file.

- `DegreeStats.jl` - This file defines the interface and different instantiations of degree statistics which are used in the estimation procedure.

- `ExactSizeCalculator.jl` - This file contains an implementation of exact subgraph counting, used to obtain the actual cardinality results for a given `QueryGraph` and `DataGraph`

- `PropertyGraph.jl` - This file defines the `QueryGraph` and `DataGraph` structs as well as functions used to edit their labels and edges.

- `QuasiStableCardinalityEstimator.jl` - This file includes an implementation of cardinality estimation using a `ColorSummary`.

- `UpdateSummary.jl` - This file includes an implementation of a proposed methodology for updating a `ColorSummary` without referencing the `DataGraph` by adjusting the `DegreeStatistics` corresponding to the given input.

- `Utils.jl` - This file includes implementations of utility methods used to find spanning trees and node orders (to decide how to traverse through a `QueryGraph`)