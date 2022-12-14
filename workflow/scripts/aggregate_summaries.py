import sys
import pandas as pd
from functools import reduce

sys.stderr = open(snakemake.log[0], "w")

# Read the summary csv files of each of the isolate processed and transposed the long table into a wide table
loaded_df = [pd.read_csv(csv, index_col=0).T for csv in snakemake.input]

# Concatenate all the genomes
aggregate_df = pd.concat(loaded_df)

# Make sure the isolate column does not remain in the index
aggregate_df = aggregate_df.rename_axis('isolate').reset_index()

# Write the aggregated table to file
aggregate_df.to_csv(snakemake.output[0], index=False)
