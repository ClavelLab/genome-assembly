rule remove_small_contigs:
    input:
        "results/assembly/{isolate}/contigs.fasta",
    output:
        "results/quality_check/{isolate}.genome.fa",
    log:
        "logs/quality_check/{isolate}_remove_small_contigs.log",
    conda:
        "../envs/seqkit.yaml"
    threads: config["threads"]
    shell:
        """
        seqkit seq --remove-gaps --min-len 1000 --threads {threads} {input} 1> {output} 2> {log}
        """


rule checkM_for_quality:
    input:
        expand(
            "results/{isolate}.{genomic_type}.fa",
            isolate=samples["isolate"],
            genomic_type=["genome", "plasmids"],
        ),
    output:
        "results/quality_check/checkm.tsv",
    log:
        "logs/quality_check/checkm.log",
    conda:
        "../envs/checkm.yaml"
    threads: config["threads"]
    shell:
        """
        checkm lineage_wf --reduced_tree \
        --extension fa \
        --tab_table --file {output} \
        --threads {threads} \
        results/ \
        results/quality_check/checkm/ &> {log}
        """
