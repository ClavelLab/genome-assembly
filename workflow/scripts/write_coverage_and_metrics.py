import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Initialize a dictionary of metrics
metrics = {}

# Read the seqkit stats table for raw sequences R1 (and possibly R2)
raw_df = pd.read_table(snakemake.input['raw'], sep="\t")

# Gather metrics on raw sequences
metrics['sequence_count'] = raw_df.at[0, 'num_seqs']
metrics['basepairs_count'] = raw_df.loc[:, 'sum_len'].sum() # if paired end data
metrics['average_length'] = raw_df.at[0, 'avg_len']

# Read the seqkit stats table for trimmed/quality sequences R1 (and possibly R2)
trimmed_df = pd.read_table(snakemake.input['trimmed'], sep="\t")

# Gather metrics on trimmed/quality filtered sequences
metrics['sequence_count_qual'] = trimmed_df.at[0, 'num_seqs']
metrics['basepairs_count_qual'] = trimmed_df.loc[:, 'sum_len'].sum() # if paired end data

# Read the seqkit stats table for the assembly
assembly_df = pd.read_table(snakemake.input['assembly'], sep="\t")

# Compute the coverage as the basepairs count of raw sequence over the basepair count of the assembly
metrics['coverage'] = round(metrics['basepairs_count'] / assembly_df.at[0, 'sum_len'])

# Export the metrics dictionary as a dataframe
metrics_df = pd.DataFrame(metrics, index=[snakemake.wildcards.isolate])
metrics_df.to_csv(snakemake.output[0])

