rule assemble_after_plasmid:
    input:
        "results/trimmed/{isolate}.1.phix.noplasmid.fastq",
        "results/trimmed/pe/{isolate}.2.phix.noplasmid.fastq",
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
