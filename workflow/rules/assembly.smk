def assess_plasmid_reconstruction_success(wildcards):
    # If initial plasmid assembly was successful
    #  provide the reads w/o plasmids for genome assembly
    # Else, re-use the phiX-removed reads for assembly
    plasmid_graph = checkpoints.plasmid_reconstruction.get(**wildcards).output[0]
    if os.path.getsize(plasmid_graph) > 0:
        extracted_graph = checkpoints.plasmid_extraction.get(**wildcards).output[0]
        if os.path.getsize(extracted_graph) > 0:
            reads = [
                "results/trimmed/{isolate}.1.phix.noplasmid.fastq",
                "results/trimmed/{isolate}.2.phix.noplasmid.fastq",
            ]
            reads.append(plasmid_graph)
            reads.append(extracted_graph)
        else:
            reads = [
                "results/trimmed/{isolate}.1.phix.fastq",
                "results/trimmed/{isolate}.2.phix.fastq",
            ]
    else:
        reads = [
            "results/trimmed/{isolate}.1.phix.fastq",
            "results/trimmed/{isolate}.2.phix.fastq",
        ]
    return reads


checkpoint assemble_after_plasmid:
    input:
        assess_plasmid_reconstruction_success,
        "results/plasmid_reconstruction/{isolate}.done",
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
