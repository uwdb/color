# Experiments

## Utility Files
The following files provide convenience when running different experiments:

- `build_color_summaries.jl` - Implements a `build_experiments` function which builds summaries for all `DataGraphs` included in the specified `ExperimentParams`. These are stored as `.csv` files.

- `convert_to_gcare.jl` - Implements a `convert_dataset` function which takes a `.graph` file (following the Subgraph-Matching format) and converts it into a `.txt` file (following the G-Care format).

- `Experiments.jl` - Provides the overall Experiments harness of necessary packages and files. Also defines the timeout constant.

- `gcare_summary.py` - Code which converts G-Care's custom output format into a `DataFrame` easily parsed in Julia.

- `get_true_cardinalities.jl` - Implements a `calculate_true_cardinalities` function obtains the true cardinality for a given dataset from a cardinality file, or calculates the exact cardinality if the file does not exist. Also has a function to verify true cardinalities.

- `graph_results.jl` - Provides different ways to graph results using the output `DataFrame` from an experiment.

- `load_datasets.jl` - Implements a `load_dataset` function which creates the appropriate `DataGraph` for the given dataset.

- `load_querysets.jl` - Implements a `load_querysets` function which creates a dictionary mapping each input dataset to its corresponding list of `QueryGraph`s.

- `run_estimators.jl` - Implements a `run_estimation_experiments` function which for each given `ExperimentParams` creates a `.csv` file describing the results after running the experiment.

- `run_submitted_experiments.sh` - This is a bash script used to run the experiments described in the submitted paper and save the corresponidng figures.

- `utils.jl` - Defines enums for different datasets as well as the `ExperimentParams` used to describe parameter settings for different experiments.

Generally, to run an experiment, one simply needs to create `ExperimentParams` describing the parameters for the experiment, call `build_experiments` to build the summaries used in the experiment, call `run_estimation_experiments` to obtain the results for the estimation, then choose an appropriate function from `graph_results.jl` to display the results. Multiple examples of this process are included in `Experiments/Scripts`.

## Scripts
There are a variety of Julia scripts provided in the `Experiments/Scripts` folder which perform different experiments then save figures presenting the results to the `Experiments/Results/Figures` folder.

For example, to find how using different maximum cycle sizes in the summary affects the cardinality estimation, the `max_cycle_size.jl` script can be called from the main directory:
```
$ julia Experiments/Scripts/max_cycle_size.jl
```

## Results
The `Experiments/Results` folder contains a folder for figures, but also any saved `ColorSummary` or cardinality estimation result from using the experiment utility files will be saved here.

### Saved `ColorSummary` Information
After calling `build_experiments`, for each `ExperimentParams`, a `.csv` file storing information about the `ColorSummary`-building process is stored in `Experiments/Results`. After the first line with the header, each row of the file stores the dataset, partitioner, number of colors, build phase, time spent in that build phase, and the overall memory footprint. These values are used in figures evaluating the performance of the framework.

Meanwhile, the actual serialized summaries are stored in `Experiments/SerializedSummaries` to be used in `run_estimation_experiments`.


### Saved Cardinality Estimation Format
After calling `run_estimation_experiments`, an estimation result `.csv` file will be saved in `Experiments/Results`. After the first line with the header, each row of the file stores the cardinality estimate, true cardinality, estimation time, query type (if included in the given query file), file path to the query file, boolean representing if the estimation failed to compute, and the path width of the query. These values are used in figures evaluating the accuracy of the framework.
