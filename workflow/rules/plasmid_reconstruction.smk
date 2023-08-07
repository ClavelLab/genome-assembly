rule dummy_plasmid:
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    shell:
        """
        cat /dev/null > {output}
        """

rule compute_lengths_plasmids:
    output:
        "results/plasmid_reconstruction/{isolate}/plasmids_lengths.tsv",
    log:
        "logs/plasmid_reconstruction/{isolate}_compute_lengths_plasmids.log",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        (seqkit fx2tab --name --length --threads {threads} {input} 1> {output} 2> {log} ) || (cat /dev/null > {output})
        """
