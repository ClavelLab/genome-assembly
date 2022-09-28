import sys
import pandas as pd
from functools import reduce
from datetime import date
import subprocess

sys.stderr = open(snakemake.log[0], "w")

# Read the csv files of basepairs metrics and SSU/LSU extraction
dict_csv = { key: pd.read_csv(snakemake.input[key], index_col=0) for key in ['metrics', 'ssu_lsu', 'trnas_5s'] }
# Read the tsv files of sample table provided in the configuration file and outputs from genome quality assessment tools
dict_table = { key: pd.read_table(snakemake.input[key], sep='\t') for key in ['samples', 'checkm', 'quast'] }
# Read the checksums of the raw fastq files and the genome fasta file
dict_md5 = { key: pd.read_table(snakemake.input[key], sep='\s+', names=['md5', 'file']) for key in ['fastq_md5', 'genome_md5'] }


# Keep the csv as a list
list_csv = list(dict_csv.values())
# Append the samples table and add to the csv list of genome metrics, SSU/LSU and tRNAs
dict_table['samples'].set_index('isolate', inplace=True)
list_csv.append(dict_table['samples'])

# Merge the table based on the index that is the isolate name
merged = reduce(lambda  left,right: pd.merge(left,right, left_index=True, right_index=True), list_csv)

# Add quality metrics from the dictionaries containing the tables
merged['compl_score'] = float(dict_table['checkm'].loc[ dict_table['checkm']['Bin Id']==snakemake.wildcards.isolate+'.genome', 'Completeness'])
merged['compl_software'] = 'checkm'
merged['contam_score'] = float(dict_table['checkm'].loc[ dict_table['checkm']['Bin Id']==snakemake.wildcards.isolate+'.genome', 'Contamination'])
merged['contam_software'] = 'checkm'

# Set the index for the QUAST data
dict_table['quast'].set_index('Assembly', inplace=True)

# Add assembly metrics
merged['N50'] = dict_table['quast'].loc[snakemake.wildcards.isolate+'.final', 'N50']
merged['genome_length'] = dict_table['quast'].loc[snakemake.wildcards.isolate+'.final', 'Total length']
merged['number_contig'] = dict_table['quast'].loc[snakemake.wildcards.isolate+'.final', '# contigs (>= 0 bp)']
merged['number_contig_below_1kb'] = dict_table['quast'].loc[snakemake.wildcards.isolate+'.raw', '# contigs (>= 0 bp)'] - dict_table['quast'].loc[snakemake.wildcards.isolate+'.raw', '# contigs (>= 1000 bp)']
merged['max_contig_length'] =  dict_table['quast'].loc[snakemake.wildcards.isolate+'.final','Largest contig']

# Add software information
merged['SSU_recover_software'] = merged['LSU_recover_software'] = 'metaxa2'
merged['trna_ext_software'] = 'tRNAscan-SE'
merged['assembly_software'] = 'SPAdes'

# Add the location and checksums of the genome
merged['genome_file'] = dict_md5['genome_md5'].at[0, 'file']
merged['genome_file_md5'] = dict_md5['genome_md5'].at[0, 'md5']

# Add information about the raw sequences: location, checksums, metrics
merged['forward_file_md5'] = dict_md5['fastq_md5'].at[0, 'md5']
merged['reverse_file_md5'] = dict_md5['fastq_md5'].at[1, 'md5']


# Rename columns
merged.rename(columns={'SSU_bp':'16S_SSU_rRNA_length', 'LSU_bp':'23S_LSU_rRNA_length', 'tRNAs':'trnas'}, inplace=True)


# Compute criteria of the MIMAG and SeqCode
hq_criteria = {
    'is_compl_grtr_90': merged.at[snakemake.wildcards.isolate, 'compl_score'] > 90,
    'is_contam_less_5': merged.at[snakemake.wildcards.isolate, 'contam_score'] < 5,
    'is_coverage_grtr_10': merged.at[snakemake.wildcards.isolate, 'coverage'] > 10,
    'are_contigs_less_100': merged.at[snakemake.wildcards.isolate, 'number_contig'] < 100,
    'is_N50_grtr_25kb': merged.at[snakemake.wildcards.isolate, 'N50'] > 25000,
    'is_max_contig_grtr_100kb': merged.at[snakemake.wildcards.isolate, 'max_contig_length'] > 100000,
    'is_trnas_grtr_18': merged.at[snakemake.wildcards.isolate, 'trnas'] > 18,
    'is_SSU_grtr_0': merged.at[snakemake.wildcards.isolate, '16S_SSU_rRNA_length'] != "0",
    'is_LSU_grtr_0': merged.at[snakemake.wildcards.isolate, '23S_LSU_rRNA_length'] != "0",
    'is_5S_grtr_0': merged.at[snakemake.wildcards.isolate, '5S_rRNA_length'] != "0",
}

# Assess which criteria of the MIMAG and SeqCode did not pass
to_check = [x for x in hq_criteria.keys() if not hq_criteria[x]]
criteria_false = ', '.join(to_check)
# Indicate the assembly quality if compliance to the MIMAG and SeqCode criteria
merged['assembly_qual'] = 'High-quality draft' if all(hq_criteria.values()) else 'Manual review: ' + criteria_false

# Add the flags of the criteria to the main table
merged = pd.merge(merged, pd.DataFrame(hq_criteria, index = [snakemake.wildcards.isolate]),
                  left_index=True, right_index=True)

# Add the plasmids lengths if a file was provided
if not snakemake.input['plasmids']:
    try:
        # Read the report
        plasmids = pd.read_table(snakemake.input['plasmids'], names = ['file', 'length'])
        # Extract the lengths and sort in descending order
        plasmids = plasmids.loc[:,'length'].sort_values(ascending = False).tolist()
        # Format the lengths as string and concatenate. No concatenation if only one element
        merged['plasmid_length'] = ';'.join([ str(x) for x in plasmids ])
    except pd.errors.EmptyDataError: # if the report is empty
        merged['plasmid_length'] = '0'
else:
    # Add a length of 0bp if no file was provided
    merged['plasmid_length'] = '0'


# Add the adapters file used
merged['adapters_file'] = snakemake.params.adapters

# Reorder columns
genome_csv = merged.reindex(columns=['genome_file', 'genome_file_md5',
                'assembly_qual','genome_length',
                'number_contig', 'N50',
                'number_contig_below_1kb', 'max_contig_length',
                'coverage', 'assembly_software',
                'compl_score', 'compl_software',
                'contam_score', 'contam_software',
                '16S_SSU_rRNA_length', 'SSU_recover_software',
                '23S_LSU_rRNA_length', 'LSU_recover_software',
                'trnas', 'trna_ext_software', '5S_rRNA_length',
                'forward_file', 'forward_file_md5',
                'reverse_file', 'reverse_file_md5',
                'sequence_count', 'basepairs_count', 'average_length',
                'sequence_count_qual', 'basepairs_count_qual', 'adapters_file', 'plasmid_length',
                'is_compl_grtr_90', 'is_contam_less_5',
                'is_coverage_grtr_10', 'are_contigs_less_100',
                'is_N50_grtr_25kb','is_max_contig_grtr_100kb',
                'is_trnas_grtr_18',
                'is_SSU_grtr_0', 'is_LSU_grtr_0', 'is_5S_grtr_0'])

# Add the date when the assembly was done
genome_csv['assembly_date'] = date.today()
# Get the version of the workflow (and the last commit if no tags are present) and add to the table
try:
    get_tag = subprocess.check_output("git describe --always", shell=True, text=True)
    genome_csv['workflow_version'] = get_tag.strip()
except subprocess.CalledProcessError:
    genome_csv['workflow_version'] = "NA"
    print("Either git is not installed or the workflow is not run in the git repository")

# Write the table to file
genome_csv.T.to_csv(snakemake.output[0])
