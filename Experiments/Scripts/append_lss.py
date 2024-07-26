# The LSS results were established using a different code base. This file takes the resulting csv files and appends the information
# to the end of an overall comparison results file to be used for plotting data in figures used in the paper.
import csv
import pandas as pd

lss_datasets = ['aids', 'dblp', 'eu2005', 'human', 'lubm80', 'yeast', 'youtube']

build_inference_filename = 'Experiments/Results/LSS/result/build_and_inference_ST.csv'
comparison_filename = 'Experiments/comparison_results.csv'

estimator = 'lss'
with open(comparison_filename, 'a') as comparison_file_obj:
    writer_obj = csv.writer(comparison_file_obj)
    for dataset in lss_datasets:
        runtime = 0
        results_filename = 'Experiments/Results/LSS/result/' + dataset + '/' + dataset + '_NNGINConcat_freq_80_cv.csv'
        with open(build_inference_filename) as build_file_obj:
            build_reader_obj = csv.reader(build_file_obj)
            for row in build_reader_obj:
                if row[0] == dataset:
                    runtime = row[3]
            build_file_obj.close()

        with open(results_filename) as results_file_obj:
            reader_obj = csv.reader(results_file_obj)
            for row in reader_obj:
                if row[2] == 'error':
                    continue
                query = 'query' + row[1]
                value = row[2]
                # now write estimator,dataset,query,value,runtime to a new row in the comparison_results.csv
                new_row = [estimator, dataset, query, value, runtime]
                writer_obj.writerow(new_row)
            results_file_obj.close()
    comparison_file_obj.close()

df = pd.read_csv('Experiments/comparison_results.csv')
df.to_parquet('Experiments/comparison_results.parquet')