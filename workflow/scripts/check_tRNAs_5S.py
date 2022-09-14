import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Read Bakta tabular annotations
annotations = pd.read_table(snakemake.input[0], sep="\t", header=2)
# Casting the column as string to fix a AttributeError when no gene names are found
#   Can only use .str accessor with string values
annotations['Gene']=annotations['Gene'].astype(str)
# Keep a subset with the tRNAs in the predicted gene name and replace NaN by False
annotations_trnas = annotations[annotations['Gene'].str.contains('_trna', na=False)]

# Dictionary counts of the number of occurrences of each provided tRNA
trnas = {x: annotations_trnas['Gene'].str.count(x+'_trna').sum() for x in snakemake.params.tRNAs}

# How many tRNAs were detected?
how_many = sum([trnas[x]>0 for x in trnas.keys()])

# Count how many 5S rRNA genes are predicted
five_s = int(annotations['Gene'].str.count('5S_rrna').sum())

# Export the detailed counts of tRNAs and 5S as a dataframe
detailled_df = pd.merge(
        pd.DataFrame(trnas, index=[snakemake.wildcards.isolate]),
        pd.DataFrame({'5S_rRNA': five_s}, index=[snakemake.wildcards.isolate]),
        left_index=True, right_index=True)
detailled_df.to_csv(snakemake.output.details)


# Export the summary count dictionary of tRNAs and 5S as a dataframe
summary_df = pd.merge(
        pd.DataFrame({'tRNAs': how_many}, index=[snakemake.wildcards.isolate]),
        pd.DataFrame({'5S_rRNA': five_s}, index=[snakemake.wildcards.isolate]),
        left_index=True, right_index=True)
summary_df.to_csv(snakemake.output.summary)

