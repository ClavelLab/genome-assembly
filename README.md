# Snakemake workflow: genome-assembly-clavel

A Snakemake workflow for assembling bacterial genomes according to the standard operating procedure in the Clavel Lab and assessing their quality in accordance to MIMAG and SeqCode criteria.

## Usage

First, the different conda environments can be installed using the following command:

    snakemake --use-conda --conda-create-envs-only

Then, the workflow can be used with the default configuration file and 9 threads as follows:

    snakemake -c 9 --use-conda

Or using a custom one:

    snakemake --configfile config/er2.yaml -c 9 --use-conda


## Input/Output

The [documentation](config/README.md) provides more details on the configuration of the workflow, but in summary:

**Input**:

* A tab-separated table indicating the isolate name and the path to paired-end FASTQ files
* A configuration file

**Output**:

* The genomes of the isolates in a gzipped FASTA file format in `results/genome/isolate.genome.fa.gz`
* The plasmids sequences of the isolates in a gzipped FASTA file format in `results/plasmids/isolate.plasmids.fa.gz`
* A comma-separated table (in `results`) summarising the metadata on the generated genomes with for instance:
    - paths to the previous files and their md5sums
    - genome quality metrics and the software used to calculate them (e.g., completion, contamination, tRNAs, 16S, 23S, 5S)
    - boolean flags indicating which genome quality tests passed or failed

## Changelog

* v3.2: Make sure to use at least the version `v1.5` of bakta
* v3.1: Reduce disk space used by flagging fast-to-generate output files as temporary (e.g., trimmed sequences, BAM file)
* v3.0: Improved README (especially for the configuration in `config`) and fixed bugs. Failed quality criteria are correctly displayed along the aggregated table, plasmids are produced only when detected and a length of 0 bp is outputted when LSU/SSU are not detected.
* v2.0: Extended the assembly quality assessment to evaluate the compliance to the MIMAG and SeqCode criteria and output an aggregated summary table
* v1.1: Added a quality control by CheckM
* v1.0: Initial conversion of the Standard Operating Procedure used in the lab to a Snakemake workflow

