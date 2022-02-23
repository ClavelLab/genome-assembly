rule plasmid_reconstruction:
    input:
        "results/trimmed/{isolate}.1.phix.fastq",
        "results/trimmed/pe/{isolate}.2.phix.fastq",
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    log:
        "logs/plasmid_reconstruction/{isolate}_plasmid_reconstruction.log",
    conda:
        "../envs/spades.yaml"
    threads: config["threads"]
    shell:
        """
        plasmidspades.py -t {threads} --careful \
            -1 {input[0]} -2 {input[1]} \
            -o results/plasmid_reconstruction/{wildcards.isolate} &> {log}
        """


rule extract_plasmid_assembly_graph:
    input:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
    log:
        "logs/plasmid_reconstruction/{isolate}_assembly_extraction.log",
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        make_fasta_from_fastg.py -g {input} -o {output} &> {log}
        """


rule index_plasmid_assembly_graph:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
    output:
        idx=multiext(
            "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/plasmid_reconstruction/{isolate}_index_assembly_graph.log",
    wrapper:
        "v1.2.0/bio/bwa/index"


rule align_plasmid_assembly_graph:
    input:
        read=[
            "results/trimmed/{isolate}.1.phix.fastq",
            "results/trimmed/pe/{isolate}.2.phix.fastq",
        ],
        idx=multiext(
            "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe.bam",
    params:
        # No sorting. But defaults include header and compressed BAM
        sort == "none",
    log:
        "logs/plasmid_reconstruction/{isolate}_align_assembly_graph.log",
    wrapper:
        "v1.2.0/bio/bwa/mem"
