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


rule extract_assembly_graph:
    input:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
    log:
        "logs/plasmid_reconstruction/{isolate}_assembly_extraction.log",
    conda:
        "../envs/recycler.yaml"
    threads: config["threads"]
    shell:
        """
        make_fasta_from_fastg.py -g {input} -o {output} &> {log}
        """
