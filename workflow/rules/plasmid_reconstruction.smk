rule plasmid_reconstruction:
    input:
        "results/trimmed/{isolate}.1.phix.fastq",
        "results/trimmed/pe/{isolate}.2.phix.fastq",
    output:
        directory("results/plasmid_reconstruction/{isolate}"),
    log:
        "logs/plasmid_reconstruction/{isolate}.log",
    conda:
        "../envs/spades.yaml"
    threads: config["threads"]
    shell:
        """
        plasmidspades.py -t {threads} --careful \
            -1 {input[0]} -2 {input[1]} \
            -o {output} &> {log}
        """
