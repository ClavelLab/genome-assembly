# Configuration

## Isolates table

The configuration file `config/config.yaml` should indicate the path to the tabular file describing the isolates and the raw sequences to be processed.
It is specified as follows: 

    samples: config/isolates.tsv

The raw sequence files can be freely named and will be standardized later. In the example below, both names are accepted.

    isolate	forward_file	reverse_file
    Favorite_isolate	/home/rickastley/raw-sequences/Favorite_isolate_R1.fastq.gz	/home/rickastley/raw-sequences/Favorite_isolate_R1.fastq.gz	
    Slowgrower_isolate	/home/rickastley/raw-sequences/ComplexName-of-the-isolate_R1.fastq.gz	/home/rickastley/raw-sequences/ComplexName-of-the-isolate_R2.fastq.gz	


## Computing settings

The maximum number of threads to be used for parallelisation should be indicated as follow (and as well during the snakemake invocation with the `-c 9` flag):

    threads: 9

The CheckM software runs a full tree by default that uses ~40Gb of RAM. If you don't have such resources, use a reduced tree (~14Gb of RAM) by setting the following flag to `true`:

    reduced_tree: false


## Databases and external resources

Place in the `resources` folder the FASTA files describing the adapters and the phiX sequences and indicate their path with these two flags:

    adapters: resources/NexteraPE-PE.fa
    phix: resources/phix.fasta

Specify where should the databases of bakta and MDMcleaner should be found or if not where they should be stored using these two flags:

    bakta_db: /data/bakta_db
    mdmcleaner_db: /data/mdmcleaner_db


