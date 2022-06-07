import sys
import re
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Initialize a dictionary of subunit length
subunit_length = {}

# Read the first line of the FASTA sequence to extract the length of the SSU (16S) and the LSU (23S)
for subunit in ['SSU', 'LSU']:
    with open(snakemake.input[subunit], 'r') as fasta:
        # Extract the first line of the FASTA
        header = fasta.readline()
    # Find length of the sequence between brackets using regex
    matches = re.findall('.*\((\d+) bp\)', header)
    # Extract and format the length
    subunit_length[subunit+'_bp'] = int(matches[0])

# Export the length dictionary as a dataframe
subunit_length_df = pd.DataFrame(subunit_length, index=[snakemake.wildcards.isolate])
subunit_length_df.to_csv(snakemake.output[0])

