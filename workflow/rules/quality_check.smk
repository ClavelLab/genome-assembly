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
        results/quality_check/checkm/
        """
