# Snakemake workflow: genome-assembly-clavel

A Snakemake workflow for assembling bacterial genomes according to the standard operating procedure in the Clavel Lab and assessing their quality in accordance to MIMAG and SeqCode criteria.

## Usage

First, the different conda environments can be installed using the following command:

    snakemake --use-conda --conda-create-envs-only

Then, the workflow can be used with the default configuration file and 9 threads as follows:

    snakemake -c 9 --use-conda

Or using a custom one:

    snakemake --configfile config/er2.yaml -c 9 --use-conda

*Note*: The default [rerun behavior](https://github.com/snakemake/snakemake/issues/1694) of Snakemake will rerun the workflow for any possible changes in input, code, parameters, modification time or software environment. However, this implies that the workflow downstream of the plasmid extraction will be rerun *everytime* because of the presence of [checkpoints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution). In a situation where the user knows that some steps do not need to be run again (e.g., plasmid extraction, contigs curation), the flag `--rerun-triggers {mtime,params,software-env,code}` (instead of the default: `--rerun-triggers {mtime,params,input,software-env,code}`) can be added at the user's responsibility.

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

*Note*: The count of tRNAs indicated in the summary **does not** correspond to the total tRNAs found **but** to the detected number of unique *essential* tRNAs among the following: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr and Val.

## Changelog

* v5.1: Fix issue in plasmid length summary. Remove bakta compliant mode.
* v5.0: Better summary with lengths of (multiple) 16S, 23S, 5S and plasmids. Set the minimum contig size to 500bp. Remove MDMcleaner step. Fix incorrect N50 and contig number by QUAST. Note on tRNAs.
* v4.2: Handle exceptions from MDMcleaner. Fix plasmid workflow to not trigger unneeded rules.
* v4.1: Fix issue with the input function to keep the workflow flexible with the plasmid assembly
* v4.0: Fix bugs such as issue with Recycler failing with unconnected plasmid assembly graph or missing gene names in annotations. Add assembly date and workflow version to the summary table. Better documentation regarding the rerun behavior of Snakemake
* v3.2: Make sure to use at least the version `v1.5` of bakta
* v3.1: Reduce disk space used by flagging fast-to-generate output files as temporary (e.g., trimmed sequences, BAM file)
* v3.0: Improved README (especially for the configuration in `config`) and fixed bugs. Failed quality criteria are correctly displayed along the aggregated table, plasmids are produced only when detected and a length of 0 bp is outputted when LSU/SSU are not detected.
* v2.0: Extended the assembly quality assessment to evaluate the compliance to the MIMAG and SeqCode criteria and output an aggregated summary table
* v1.1: Added a quality control by CheckM
* v1.0: Initial conversion of the Standard Operating Procedure used in the lab to a Snakemake workflow

