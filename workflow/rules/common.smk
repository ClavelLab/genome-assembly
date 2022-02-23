import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("isolate", drop=False)

validate(samples, schema="../schemas/samples.schema.yaml")

# Fetch the corresponding fastq files from the given isolate
def get_fastqs(wildcards):
    return dict(zip(["r1", "r2"], samples.loc[wildcards.isolate, ["fq1", "fq2"]]))


rule rename_genome_fasta:
    input:
        "results/assembly/{isolate}/contigs.fasta",
    output:
        "results/{isolate}.genome.fa",
    shell:
        "mv {input} {output}"


rule rename_plasmid_fasta:
    input:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    output:
        "results/{isolate}.plasmids.fa",
    shell:
        "mv {input} {output}"
