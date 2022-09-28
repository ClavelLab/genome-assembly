import sys
import re
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Initialize a dictionary of subunit length
subunit_length = {}

# Read the first line of the report giving the length of all the SSU (16S) and the LSU (23S) FASTA sequences detected
for subunit in ['SSU', 'LSU']:
    try:
        # Read the report
        ssu_lsu = pd.read_table(snakemake.input[subunit], names = ['file', 'length'])
        # Extract the lengths and sort in descending order
        ssu_lsu_lengths = ssu_lsu.loc[:,'length'].sort_values(ascending = False).tolist()
        # Format the lengths as string and concatenate. No concatenation if only one element
        subunit_length[subunit+'_bp'] = ';'.join([ str(x) for x in ssu_lsu_lengths ])
    except pd.errors.EmptyDataError: # if the report is empty
        subunit_length[subunit+'_bp'] = "0"
# Export the length dictionary as a dataframe
subunit_length_df = pd.DataFrame(subunit_length, index=[snakemake.wildcards.isolate])
subunit_length_df.to_csv(snakemake.output[0])

