rule remove_small_contigs:
    input:
        "results/assembly/{isolate}/contigs.fasta",
    output:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
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
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        "results/quality_check/{isolate}/{isolate}_checkm.tsv",
    log:
        "logs/quality_check/{isolate}_checkm.log",
    params:
        flag="--reduced_tree" if config["reduced_tree"] else "",
    conda:
        "../envs/checkm.yaml"
    threads: config["threads"]
    shell:
        """
        checkm lineage_wf {params.flag} \
        --extension fa \
        --tab_table --file {output} \
        --threads {threads} \
        results/quality_check/{wildcards.isolate}/ \
        results/quality_check/{wildcards.isolate}/checkm/ &> {log}
        """


rule download_bakta_db:
    output:
        directory(config["bakta_db"]),
    log:
        "logs/quality_check/bakta_db.log",
    conda:
        "../envs/bakta.yaml"
    shell:
        """
        bakta_db download --output {output}
        """


rule bakta_for_annotation:
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        directory("results/quality_check/{isolate}/bakta/"),
    log:
        "logs/quality_check/{isolate}_bakta.log",
    params:
        db=config["bakta_db"],
    conda:
        "../envs/bakta.yaml"
    threads: config["threads"]
    shell:
        """
        bakta --db {params.db}/db/ \
        --prefix {wildcards.isolate} \
        --output {output} \
        --threads {threads} {input} &> {log}
        """


rule check_tRNAs_5S:
    input:
        "results/quality_check/{isolate}/bakta/{isolate}.tsv",
    output:
        summary="results/quality_check/{isolate}/bakta/{isolate}.tRNAs-5S.csv",
        details="results/quality_check/{isolate}/bakta/{isolate}.tRNAs-5S.details.csv",
    log:
        "logs/quality_check/{isolate}_check_tRNAs_5S.log",
    conda:
        "../envs/pandas.yaml"
    params:
        tRNAs=[
            "Ala",
            "Arg",
            "Asn",
            "Asp",
            "Cys",
            "Gln",
            "Glu",
            "Gly",
            "His",
            "Ile",
            "Leu",
            "Lys",
            "Met",
            "Phe",
            "Pro",
            "Ser",
            "Thr",
            "Trp",
            "Tyr",
            "Val",
        ],
    script:
        "../scripts/check_tRNAs_5S.py"
