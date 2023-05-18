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


def plasmids_when_needed():
    # Trigger the complete plasmid extraction workflow
    #  if initial plasmid assembly was successful
    plasmids = []
    for iso in samples["isolate"]:
        plasmid_graph = checkpoints.plasmid_reconstruction.get(isolate=iso).output[0]
        # Check first the plasmid reconstruction step in order
        # to not trigger the plasmid workflow if the extraction was not run
        if os.path.getsize(plasmid_graph) > 0:
            # If the extracted graph was done, it means the plasmid reconstruction was done
            # but checking only the extraction could trigger the workflow when not needed
            extracted_graph = checkpoints.plasmid_extraction.get(isolate=iso).output[0]
            # If file is non empty
            if os.path.getsize(extracted_graph) > 0:
                plasmids.append(iso)
    return plasmids
