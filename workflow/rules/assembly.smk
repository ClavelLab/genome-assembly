def assess_plasmid_reconstruction_success(wildcards):
    # If initial plasmid assembly was successful
    #  provide the reads w/o plasmids for genome assembly
    # Else, re-use the phiX-removed reads for assembly
    plasmid_graph = checkpoints.plasmid_reconstruction.get(**wildcards).output[0]
    extracted_graph = checkpoints.plasmid_extraction.get(**wildcards).output[0]
    with plasmid_graph.open() as f:
        # Read the first caracter of the plasmid graph
        plasmid = f.read(1)
    with extracted_graph.open() as f:
        # Read the first caracter of the extracted graph
        extracted = f.read(1)
    # If any files are empty
    if not plasmid or not extracted:
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
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
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
