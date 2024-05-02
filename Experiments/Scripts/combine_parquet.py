import os
import pandas as pd
import re
import pathlib

import pyarrow.parquet as pq

path1 = "Experiments/comparison_results.parquet" # this should represent the path to the overall results
path2 = "Experiments/alley_results.parquet" # this should represent the path to just the alleyTPI results


df1 = pq.read_table(source=path1, filters = [('Estimator', '!=', 'alleyTPI'), ('Estimator', '!=', 'alley')]).to_pandas() # filter out the old alley results
df2 = pq.read_table(source=path2).to_pandas() # collect all the new alley results

df_result = pd.concat([df1, df2]) # combine :)

df_result.to_parquet('comparison_results.parquet')
df_result.to_csv('comparison_results.csv', index=False)
