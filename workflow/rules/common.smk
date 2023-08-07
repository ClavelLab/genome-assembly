import pandas as pd
import yaml as yaml
import os
from snakemake.utils import validate
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("isolate", drop=False)

validate(samples, schema="../schemas/samples.schema.yaml")

# Fetch the corresponding fastq files from the given isolate
def get_fastqs(wildcards):
    return dict(
        zip(
            ["r1", "r2"],
            samples.loc[wildcards.isolate, ["forward_file", "reverse_file"]],
        )
    )


rule rename_genome_fasta:
    input:
        "results/quality_check/{isolate}/{isolate}.combined.fa",
    output:
        "results/genome/{isolate}.combined.fa.gz",
    shell:
        "cat {input} | gzip > {output}"


use rule rename_genome_fasta as rename_plasmid_fasta with:
    input:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    output:
        "results/plasmids/{isolate}.plasmids.fa.gz",

