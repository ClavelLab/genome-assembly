# Snakemake workflow: genome-assembly-clavel

A Snakemake workflow for assembling bacterial genomes according to the standard operating procedure in the Clavel Lab and assessing their quality in accordance to MIMAG and SeqCode criteria.

## Usage

First, the different conda environments can be installed using the following command:

    snakemake --use-conda --conda-create-envs-only

Then, the workflow can be used with the default configuration file and 9 threads as follows:

    snakemake -c 9 --use-conda

Or using a custom one:

    snakemake --configfile config/er2.yaml -c 9 --use-conda


### Notes

1. The default [rerun behavior](https://github.com/snakemake/snakemake/issues/1694) of Snakemake will rerun the workflow for any possible changes in input, code, parameters, modification time or software environment. However, this implies that the workflow downstream of the plasmid extraction will be rerun *everytime* because of the presence of [checkpoints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution). In a situation where the user knows that some steps do not need to be run again (e.g., plasmid extraction, contigs curation), the flag `--rerun-triggers {mtime,params,software-env,code}` (instead of the default: `--rerun-triggers {mtime,params,input,software-env,code}`) can be added at the user's responsibility.
2. The assembly workflow is done by default only with the paired reads, but it is now possible to use the unpaired reads (quality-filtered, adapters- and phix-removed) when needed (e.g., low quality reverse reads that bring the paired reads number too low). To enable this feature, one must change the default git branch with the following command: `git checkout with-unpaired` and run the workflow according to the Usage section.

## Input/Output

The [documentation](config/README.md) provides more details on the configuration of the workflow, but in summary:

**Input**:

* A tab-separated table indicating the isolate name and the path to paired-end FASTQ files
* A configuration file

**Output**:

* The genome of the isolate (with the plasmids sequences if any) in a gzipped FASTA file format in `results/genome/isolate.combined.fa.gz`
* The plasmids sequences of the isolates in a gzipped FASTA file format in `results/plasmids/isolate.plasmids.fa.gz`
* A comma-separated table (in `results`) summarising the metadata on the generated genomes with for example:
    - paths to the previous files and their md5sums
    - genome quality metrics and the software used to calculate them (e.g., completion, contamination, tRNAs, 16S, 23S, 5S)
    - boolean flags indicating which genome quality tests passed or failed

The column of the output table are detailed in the next section with examples.

### Output table metadata details

| Metadata                   | Definition                                                                                                                                                                                   | Example                                      |
|----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------|
| `isolate`                  | Identifier of the isolate provided in the sample table                                                                                                                                       | CLAXXH123                                    |
| `genome_md5`               | Checksum of the FASTA file of the genome                                                                                                                                                     | db2cd4eff6d40c463998851f43b15670             |
| `assembly_qual`            | Assembly quality, either “High-quality draft” if all quality flags are valid or “Manual review:” with the quality flags that were invalid.                                                   | Manual review: are_contigs_less_100          |
| `genome_length`            | Size of the genome in nucleotide based on the length of all contigs in the assembly                                                                                                          | 4839229                                      |
| `number_contig`            | Number of contigs in the assembly                                                                                                                                                            | 134                                          |
| `N50`                      | Median length of the contigs in the assembly                                                                                                                                                 | 83304                                        |
| `number_contig_below_1kb`  | Number of small contigs, under 1000 nucleotides-long                                                                                                                                         | 202                                          |
| `max_contig_length`        | Length of the longest contig                                                                                                                                                                 | 291966                                       |
| `coverage`                 | Ratio of the basepairs count of raw sequences over the basepairs count of the assembly                                                                                                       | 394                                          |
| `assembly_software`        | Name of the software used for the assembly                                                                                                                                                   | SPAdes                                       |
| `compl_score`              | Percentage of estimated marker-based completion of the assembly                                                                                                                              | 99.46                                        |
| `compl_software`           | Name of the software used for assessing the assembly completion                                                                                                                              | checkm                                       |
| `contam_score`             | Percentage of estimated marker-based contamination of the assembly                                                                                                                           | 0                                            |
| `contam_software`          | Name of the software used for assessing the assembly contamination                                                                                                                           | checkm                                       |
| `16S_SSU_rRNA_length`      | Semicolon list of lengths of detected 16S rRNA sequences. Lengths are sorted decreasing or a value of 0 when no sequences are detected.                                                      | 1393;1037;597;528                            |
| `SSU_recover_software`     | Name of the software used for detecting the 16S rRNA sequences                                                                                                                               | metaxa2                                      |
| `23S_LSU_rRNA_length`      | Semicolon list of lengths of detected 23S rRNA sequences. Lengths are sorted decreasing or a value of 0 when no sequences are detected.                                                      | 1637;795;528                                 |
| `LSU_recover_software`     | Name of the software used for detecting the 23S rRNA sequences                                                                                                                               | metaxa2                                      |
| `trnas`                    | Number of unique essential tRNAs detected in the assembly among the following twenty: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr and Val. | 20                                           |
| `trna_ext_software`        | Name of the software used for detecting tRNAs sequences                                                                                                                                      | tRNAscan-SE                                  |
| `5S_rRNA_length`           | Semicolon list of lengths of detected 5S rRNA sequences. Lengths are sorted decreasing or a value of 0 when no sequences are detected.                                                       | 110;110;106;105                              |
| `archive_file`             | Path to the archive of the genome file                                                                                                                                                       | results/genome/CLAJMH5.genome.fa.gz          |
| `archive_file_md5`         | Checksum of the archive of the genome file                                                                                                                                                   | db2cd4eff6d40c463998851f43b15670             |
| `forward_file`             | Path to the forward/R1 raw sequence file                                                                                                                                                     | /DATA/hibc_fastq/CLAJMH5_S16_R1_001.fastq.gz |
| `forward_file_md5`         | Checksum of the archive of the forward raw sequence file                                                                                                                                     | 1ff2885fe846135980b481996d8e1da6             |
| `reverse_file`             | Path to the reverse/R2 raw sequence file                                                                                                                                                     | /DATA/hibc_fastq/CLAJMH5_S16_R2_001.fastq.gz |
| `reverse_file_md5`         | Checksum of the archive of the reverse raw sequence file                                                                                                                                     | 3acd93ba96b97aac4da7756c83f2f990             |
| `sequence_count`           | Number of reads in the library (sequencing depth)                                                                                                                                            | 6314677                                      |
| `basepairs_count`          | Number of basepairs in nucleotides                                                                                                                                                           | 1907032454                                   |
| `average_length`           | The average read length estimated as basepairs_count divided by sequence_count                                                                                                               | 151                                          |
| `sequence_count_qual`      | Number of reads in the library (sequencing depth) after quality filtering                                                                                                                    | 4837051                                      |
| `basepairs_count_qual`     | Number of basepairs in nucleotides after quality filtering                                                                                                                                   | 1324608613                                   |
| `adapters_file`            | Path to the FASTA file containing the adapters sequences                                                                                                                                     | resources/NexteraPE-PE.fa                    |
| `plasmid_length`           | Semicolon list of lengths of detected plasmid sequences. Lengths are sorted decreasing or a value of 0 when no sequences are detected.                                                       | 7059                                         |
| `is_compl_grtr_90`         | Logical flag indicating whether the completion, in compl_score, is greater than 90%                                                                                                          | TRUE                                         |
| `is_contam_less_5`         | Logical flag indicating whether the contamination, contam_score,  is greater than 5%                                                                                                         | TRUE                                         |
| `is_coverage_grtr_10`      | Logical flag indicating whether the coverage is greater than 10x                                                                                                                             | TRUE                                         |
| `are_contigs_less_100`     | Logical flag indicating whether the number of contigs, in number_contig, is below 100                                                                                                        | FALSE                                        |
| `is_N50_grtr_25kb`         | Logical flag indicating whether the N50 is longer than 25 000 nucleotides                                                                                                                    | TRUE                                         |
| `is_max_contig_grtr_100kb` | Logical flag indicating whether the longest contig, in max_contig_length, is longer than 100 000 nucleotides                                                                                 | TRUE                                         |
| `is_trnas_grtr_18`         | Logical flag indicating whether the number of tRNAs is greater than 18                                                                                                                       | TRUE                                         |
| `is_SSU_grtr_0`            | Logical flag indicating whether any 16S/SSU rRNA sequence longer than 0 nucleotide is detected                                                                                               | TRUE                                         |
| `is_LSU_grtr_0`            | Logical flag indicating whether any 23S/LSU rRNA sequence longer than 0 nucleotide is detected                                                                                               | TRUE                                         |
| `is_5S_grtr_0`             | Logical flag indicating whether any 5S rRNA sequence longer than 0 nucleotide is detected                                                                                                    | TRUE                                         |
| `assembly_date`            | Date on which the assembly was generated                                                                                                                                                     | 2022-10-12                                   |
| `workflow_version`         | Git tag of the genome assembly workflow used. If not a round version like v5.2, then the commit hash past the round version is added                                                         | v5.2-1-g4dae043                              |



## Changelog

* v7.0: Concatenate genome and plasmids sequences (if any). Add plasmids proteins annotations.
* v6.1: Fix conda environments setup by updating snakemake wrappers.
* v6.0: Checksums for both the genome archive and FASTA file are computed. Improved the table. Assembly with unpaired reads is possible.
* v5.3.2: Fixes for conda environments
* v5.3: Update bakta (fix protein predictions) and minor updates ofcheckm, metaxa2 and SPades. Better documentation of the column of the output table.
* v5.2: Recovered the behavior of length of 0bp when LSU/SSU/plasmids are not detected.
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

