# Configuration

## Isolates table

The configuration file `config/config.yaml` should indicate the path to the tabular file describing the isolates and the raw sequences to be processed.
It is specified as follows: 

    samples: config/isolates.tsv

The raw sequence files can be freely named and will be standardized later. In the example below, both names are accepted.

    isolate	forward_file	reverse_file
    Favorite_isolate	/home/rickastley/raw-sequences/Favorite_isolate_R1.fastq.gz	/home/rickastley/raw-sequences/Favorite_isolate_R1.fastq.gz	
    Slowgrower_isolate	/home/rickastley/raw-sequences/ComplexName-of-the-isolate_R1.fastq.gz	/home/rickastley/raw-sequences/ComplexName-of-the-isolate_R2.fastq.gz	

### An example of semi-automatic creation of the table of isolates and FASTQ files

In the configuration directory of the workflow:

    cd config/

List all the forward (R1) and reverse (R2) FASTQ files in their directory

    ls /data/project_fastq/PROJECT*R1_001.fastq.gz > R1
    ls /data/project_fastq/PROJECT*R2_001.fastq.gz > R2

Extract the isolate name located before the first underscore

    # Example of the path /data/PROJECT_fastq/PROJECTH117_S106_R1_001.fastq.gz
    sed -e 's/^.*PROJECT/PROJECT/' -e 's/_.*//' R1 > sample

Paste together the columns and separate with tabulations (default)

    paste -- sample R1 R2 > project

Add the necessary header to generate the required table as `project.tsv`

    cat <(echo -e "isolate\tforward_file\treverse_file" ) project > project.tsv
    rm R1 R2 sample project

*Note*: Don't forget to update the configuration file with the correct filename!

## Computing settings

The maximum number of threads to be used for parallelisation should be indicated as follow (and as well during the snakemake invocation with the `-c 9` flag):

    threads: 9

The CheckM software runs a full tree by default that uses ~40Gb of RAM. If you don't have such resources, use a reduced tree (~14Gb of RAM) by setting the following flag to `true`:

    reduced_tree: false


## Databases and external resources

Place in the `resources` folder the FASTA files describing the adapters and the phiX sequences and indicate their path with these two flags:

    adapters: resources/NexteraPE-PE.fa
    phix: resources/phix.fasta

Specify where should the database of bakta should be found or if not where they should be stored using this flag:

    bakta_db: /data/bakta_db

## Assembly parameters

Set the minimum contig size to be kept in the assembly in basepair:

    min_contig_length: 500
