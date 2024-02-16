# Code to convert GCare's custom output format into a DataFrame 

import os
import pandas as pd
import re
import pathlib

# Assuming you are in the directory where the files are located
directory = '.'  # Change this to your specific directory

# Prepare an empty list to store the records
records = []

# Regex to parse filename and extract dataset, algorithm, etc.
filename_pattern = re.compile(r'(\w+)_([a-z]+)_p(\d\.\d+)_s(\d)_query_result\.txt')

# Scan the directory for .txt files and ignore .err files
for filename in os.listdir(directory):
    if filename.endswith('.txt') and not filename.endswith('.err'):
        # Extract information from the filename
        match = filename_pattern.match(filename)
        if match:
            dataset, algorithm, probability, seed = match.groups()

            # Construct the file path
            file_path = os.path.join(directory, filename)

            # Open and read the contents of the file
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for i in range(0, len(lines), 2):  # Process files in pairs of lines
                    query_line = lines[i].strip()
                    if query_line[0] == '#':
                        continue
                    result_line = lines[i+1].strip()

                    # Extract the query from the query line
                    metadata = query_line.split()
                    query_path = pathlib.Path(metadata[0])
                    query = '/'.join(query_path.parts[-2:])
                    runtime = float(metadata[2])

                    # Extract runtime and result(estimation) values
                    values = result_line.split()
                    estimate = float(values[1])  # All subsequent values are estimates

                    # Record multiple estimates for the same query and dataset
                    records.append({
                        'Estimator': algorithm,
                        'Dataset': dataset,
                        'Query': query,
                        'Value': estimate,
                        'Runtime': runtime
                    })

# Create the dataframe
df = pd.DataFrame(records)

# The dataframe with schema: Estimate(Estimator, Dataset, Query, Value, Runtime)
print(df.head())  # Print the first few rows to check

# Save the dataframe to a CSV file if needed
df.to_csv('experiment_results.csv', index=False)
df.to_parquet('experiment_results.parquet')
