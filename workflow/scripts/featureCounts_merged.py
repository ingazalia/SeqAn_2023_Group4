import os
import pandas as pd

samples_small = pd.read_csv(str(snakemake.input["sheet"]), index_col=False,sep='\t')

file = pd.read_csv("results/featureCount/"+str(samples_small['sample'][1]) + "_counts.txt", index_col='Geneid',sep='\t',header=1) #TODO: added str

merged = pd.concat([file.iloc[:,-1]], axis=1)

for names in samples_small['sample'][1:]:
    file2 = pd.read_csv("results/featureCount/"+str(names) + "_counts.txt", index_col='Geneid',sep='\t',header=1) #TODO: added str
    merged = pd.concat([merged,file2.iloc[:,-1]], axis=1)

row_strings = (file[(file.columns)[0:-2]]).apply(lambda row: '_'.join(row.astype(str)), axis=1)

merged = pd.concat([row_strings,merged], axis=1)

# new column names
column_names = (samples_small["sample"].values).tolist()
column_names.insert(0,'Gene')

# Specify the output file path
output_file = str(snakemake.output)

# Iterate over the rows of the DataFrame
with open(output_file, 'w') as file:
    file.write('\t'.join(str(value) for value in column_names) + '\n')
    for _, row in merged.iterrows():
        file.write('\t'.join(str(value) for value in row) + '\n')