def assess_plasmid_reconstruction_success(wildcards):
    # If plasmid assembly was successful
    #  trigger the complete plasmid extraction workflow
    # Else, re-use the phiX-removed reads for assembly
    plasmid_graph = checkpoints.plasmid_reconstruction.get(**wildcards).output[0]
    with plasmid_graph.open() as f:
        # Read the first caracter of the graph
        plasmid = f.read(1)
    # If file is empty
    if not plasmid:
        message("No plasmid detected in the raw sequences. Moving to assembly")
        reads = [
            "results/trimmed/{isolate}.1.phix.fastq",
            "results/trimmed/{isolate}.2.phix.fastq",
        ]
    else:
        reads = [
            "results/trimmed/{isolate}.1.phix.noplasmid.fastq",
            "results/trimmed/{isolate}.2.phix.noplasmid.fastq",
        ]
    # Indicate the adequate reads to use
    return reads


rule assemble_after_plasmid:
    input:
        assess_plasmid_reconstruction_success,
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    output:
        "results/assembly/{isolate}/contigs.fasta",
    log:
        "logs/assembly/{isolate}.log",
    conda:
        "../envs/spades.yaml"
    threads: config["threads"]
    shell:
        """
        spades.py -t {threads} --careful \
            -1 {input[0]} -2 {input[1]} \
            -o results/assembly/{wildcards.isolate} &> {log}
        """
