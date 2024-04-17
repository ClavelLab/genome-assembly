import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Read Bakta tabular annotations
# Header of TSV files changes between v1.6.3 and v1.9.3
# https://git.rwth-aachen.de/clavellab/genome-assembly/-/issues/14
# So more flexible input options.
annotations = pd.read_table(snakemake.input[0], sep="\t", header=0, comment = "#",
        names = ["Sequence Id","Start","Stop","Strand","Locus Tag","Gene","Product","DbXrefs"])
# Casting the column as string to fix a AttributeError when no gene names are found
#   Can only use .str accessor with string values
annotations['Gene']=annotations['Gene'].astype(str)
# Keep a subset with the tRNAs in the predicted gene name and replace NaN by False
annotations_trnas = annotations[annotations['Gene'].str.contains('_trna', na=False)]

# Dictionary counts of the number of occurrences of each provided tRNA
trnas = {x: annotations_trnas['Gene'].str.count(x+'_trna').sum() for x in snakemake.params.tRNAs}

# How many tRNAs were detected?
how_many = sum([trnas[x]>0 for x in trnas.keys()])

# Subset the annotations with the prediceted 5S rRNA genes
five_s = annotations[annotations.Gene == '5S_rrna']
# Compute the length
five_s_lengths = five_s.Stop - five_s.Start
# Extract the lengths and sort in descending order
five_s_lengths = five_s_lengths.sort_values(ascending = False).tolist()
# Format the lengths as strings and concatenate. No concatenation if only one element
if five_s_lengths:
    five_s_lengths = ';'.join([ str(x) for x in five_s_lengths ])
else:
    five_s_lengths = '0'



# Export the detailed counts of tRNAs as a dataframe
detailled_df = pd.DataFrame(trnas, index=[snakemake.wildcards.isolate])
detailled_df.to_csv(snakemake.output.details)


# Export the summary count dictionary of tRNAs and 5S as a dataframe
summary_df = pd.merge(
        pd.DataFrame({'tRNAs': how_many}, index=[snakemake.wildcards.isolate]),
        pd.DataFrame({'5S_rRNA_length': five_s_lengths}, index=[snakemake.wildcards.isolate]),
        left_index=True, right_index=True)
summary_df.to_csv(snakemake.output.summary)

