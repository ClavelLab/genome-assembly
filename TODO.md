## Configuration and flexibility

- [ ] Set a flexible configuration file to be able to assemble:
    - [ ] one genome 
    - [x] multiple genomes
- [ ] Set the necessary adaptater FASTA files depending on the technology (NextSeq or MiSeq) or allow detection from filename
- [ ] The adapter files can be downloaded from the [trimmomatic repository](https://github.com/timflutre/trimmomatic/tree/master/adapters) and the PhiX genome (`NC_001422.1`) as well to enable full reproducibility

## Reporting and quality control

- [ ] Include a report rule
- [ ] Assess assembly quality with checkM. Available in bioconda (1.1.3)
- [ ] Compile assembly stats with Tom Hitch script from his github. Possibly use the `github()` directive and/or the auxiliary source files

### Alternatives to quality control

- [ ] Assess assembly quality with BUSCO (contamination from euk; more efficient?). Available in bioconda (5.2.2)
- [ ] Assess assembly quality with Merqury 
- [ ] Assess assembly quality with GUNC that can be combined natively with checkM to provide uncluttered genomes 


## Convert the SOP steps into a Snakemake workflow

Listed in the reverse order because it is easier for Snakemake design. The subsections could serve as building separate Snakefiles to be included.

### File management

- [ ] Generate genome FASTA file only
- [ ] Generate genome FASTA with plasmids if present (consider snakemake checkpoints for evaluation of condition)

### Genome assembly

- [ ] Assemble with spades (v3.13.1). Snakemake wrapper only for metaspades. Available in bioconda (3.15.3)

### Plasmid reconstruction

- [ ] Remove plasmid contigs from reads with bbduk included in the bbmap (v38.84). Snakemake wrapper available (38.90)
- [ ] Extract plasmid sequences with recycler (v unknowm) from de novo assembly graph and alignment. Available in bioconda (v0.7)
- [ ] BAM/SAM management with samtools (v0.1.19). Snakemake wrapper available (1.10)
- [ ] Alignement of reads on the assembly graph with bwa mem (v0.7.5). Snakemake wrapper available (0.7.17)
- [ ] Indexing of the assembly graph with bwa (v0.7.5). Snakemake wrapper available (0.7.17)
- [x] Convert the assembly graph in FASTA with `make_fasta_from_fastg` from Recycler (0.62). Available in bioconda (0.7-3)
- [x] Plasmid reconstruction with plasmidspades (v3.13.1). Snakemake wrapper only for metaspades. Available in bioconda (3.15.3)

### Quality filtering 

- [x] Remove phiX sequences from reads with bbduk included in the bbmap (v38.84). Snakemake wrapper available
- [x] Remove adaptaters and filter length with trimmomatic (v0.39). Snakemake wrapper available but older (0.36) so bioconda.
 
