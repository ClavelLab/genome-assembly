## Configuration and flexibility

- [ ] Set a flexible configuration file to be able to assemble:
    - [ ] one genome 
    - [x] multiple genomes
- [ ] Set the necessary adapter FASTA files depending on the technology (NextSeq or MiSeq) or allow detection from filename
- [ ] The adapter files can be downloaded from the [trimmomatic repository](https://github.com/timflutre/trimmomatic/tree/master/adapters) and the PhiX genome (`NC_001422.1`) as well to enable full reproducibility
- [ ] SPAdes provide an alternative flag to `--careful` that is `--isolate` (introduced in 3.14.0) that [could be used](https://github.com/ablab/spades/blob/spades_3.15.4/README.md#sec3.2) for high-coverage (100x) isolate genome. Note that there is no one-size-fits-all as [always](https://github.com/ablab/spades/issues/600)
- [ ] Same for recycler with a `-i True` for isolate
- [x] **Adjust the maximum length of the kmer required by recycler based on the SPAdes output**

## Reporting and quality control

- [x] Compute the numbers of contigs below 1kb and remove with seqkit
- [x] Assess completeness and contamination with checkM, but remove the plasmid check. Available in bioconda (1.1.3)
- [x] Annotate genome with Bakta (5S extraction)
- [x] Extract the LSU (23S) and SSU (16S) with metaxa2
- [x] Assess contamination with MDMcleaner
- [x] Compute basepairs statistics and coverage with seqkit
- [x] Compute assembly statistics with QUAST
- [ ] Generate checksums with md5 hash on the gz version of the raw reads and the final genome for deposition on Coscine
- [ ] Extract above statistics to produce a standard compliant table based on the sample table
- [ ] Include a report rule


## Convert the SOP steps into a Snakemake workflow

Listed in the reverse order because it is easier for Snakemake design. The subsections could serve as building separate Snakefiles to be included.

### File management

- [ ] Generate genome FASTA file only
- [x] Generate genome FASTA with plasmids if present (consider snakemake checkpoints for evaluation of condition)

### Genome assembly

- [x] Assemble with spades (v3.13.1). Snakemake wrapper only for metaspades. Available in bioconda (3.15.3)

### Plasmid reconstruction

- [x] Remove plasmid contigs from reads with bbduk included in the bbmap (v38.84). Snakemake wrapper available (38.90)
- [x] Extract plasmid sequences with recycler (v unknowm) from de novo assembly graph and alignment. Available in bioconda (v0.7)
- [x] BAM/SAM management with samtools (v0.1.19). Snakemake wrapper available (1.10)
- [x] Alignement of reads on the assembly graph with bwa mem (v0.7.5). Snakemake wrapper available (0.7.17)
- [x] Indexing of the assembly graph with bwa (v0.7.5). Snakemake wrapper available (0.7.17)
- [x] Convert the assembly graph in FASTA with `make_fasta_from_fastg` from Recycler (0.62). Available in bioconda (0.7-3)
- [x] Plasmid reconstruction with plasmidspades (v3.13.1). Snakemake wrapper only for metaspades. Available in bioconda (3.15.3)

### Quality filtering 

- [x] Remove phiX sequences from reads with bbduk included in the bbmap (v38.84). Snakemake wrapper available
- [x] Remove adapters and filter length with trimmomatic (v0.39). Snakemake wrapper available but older (0.36) so bioconda.
 
