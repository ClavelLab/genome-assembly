import sys
import pandas as pd
import re

sys.stderr = open(snakemake.log[0], "w")

# Read Bakta tabular annotations
# Header of TSV files changes between v1.6.3 and v1.9.3
# https://git.rwth-aachen.de/clavellab/genome-assembly/-/issues/14
# So more flexible input options.
annotations = pd.read_table(snakemake.input[0], sep="\t", header=0, comment = "#",
        names = ["Sequence Id","Type","Start","Stop","Strand","Locus Tag","Gene","Product","DbXrefs"])
# Casting the column as string to fix a AttributeError when no gene names are found
#   Can only use .str accessor with string values
annotations['Product']=annotations['Product'].astype(str)
# Keep a subset with the tRNAs in the type but not a pseudo tRNA nor undefined
# replace NaN by True to not select them
annotations_trnas = annotations[ (annotations['Type'] == "tRNA" ) &
        ( ~ annotations['Product'].str.contains("pseudo", na = True) ) &
        ( ~ annotations['Product'].str.contains("Xxx", na = True) ) ]


# Gene: trnL
# Product: tRNA-Leu(taa)
# Extracted: Leu
pattern = "tRNA-([A-Za-z]+)"
extracted_trnas = annotations_trnas['Product'].apply(lambda x: re.findall(pattern,x)[0])

d_trnas_count = extracted_trnas.value_counts().to_dict()
# https://stackoverflow.com/a/5352630/21085566
# Dictionary counts of the number of occurrences of each provided tRNA
trnas = {aa: d_trnas_count.get(aa, 0) for aa in snakemake.params.tRNAs}

# How many tRNAs were detected?
how_many = sum([trnas[x]>0 for x in trnas.keys()])

# Subset the annotations with the predicted 5S rRNA genes
five_s = annotations[annotations.Gene == 'rrf']
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

