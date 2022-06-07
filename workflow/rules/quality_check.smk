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


rule check_SSU:
    # SSU: 16S (BactArch) 18S (Euk)
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        fasta=temp("SSU-{isolate}.extraction.fasta"),
        results=temp("SSU-{isolate}.extraction.results"),
    log:
        "logs/quality_check/{isolate}_check_SSU.log",
    conda:
        "../envs/metaxa2.yaml"
    threads: config["threads"]
    group:
        "SSU"
    shell:
        """
        metaxa2_x -i {input} \
            -o SSU-{wildcards.isolate} -f fasta \
            -g ssu --mode genome \
            --table F --fasta T --graphical F --summary F \
            --cpu {threads} &> {log}
        """


rule move_SSU_sequence:
    input:
        "SSU-{isolate}.extraction.fasta",
    output:
        "results/quality_check/{isolate}/{isolate}.SSU.fa",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    group:
        "SSU"
    shell:
        """
        seqkit replace -p 'NODE' -r '{wildcards.isolate} SSU NODE' --line-width 0 --threads {threads} {input} > {output}
        """


rule check_LSU:
    # LSU: 23 (BactArch) 5.8 & 28S (Euk)
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        fasta=temp("LSU-{isolate}.extraction.fasta"),
        results=temp("LSU-{isolate}.extraction.results"),
    log:
        "logs/quality_check/{isolate}_check_LSU.log",
    conda:
        "../envs/metaxa2.yaml"
    threads: config["threads"]
    group:
        "LSU"
    shell:
        """
        metaxa2_x -i {input} \
            -o LSU-{wildcards.isolate} -f fasta \
            -g lsu --mode genome \
            --table F --fasta T --graphical F --summary F \
            --cpu {threads} &> {log}
        """


rule move_LSU_sequence:
    input:
        "LSU-{isolate}.extraction.fasta",
    output:
        "results/quality_check/{isolate}/{isolate}.LSU.fa",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    group:
        "LSU"
    shell:
        """
        seqkit replace -p 'NODE' -r '{wildcards.isolate} LSU NODE' --line-width 0 --threads {threads} {input} > {output}
        """


rule aggregate_SSU_LSU_results:
    input:
        SSU="SSU-{isolate}.extraction.results",
        LSU="LSU-{isolate}.extraction.results",
    output:
        wd=directory("results/quality_check/{isolate}/metaxa"),
        summary="results/quality_check/{isolate}/metaxa/{isolate}.SSU-LSU.csv",
    shell:
        """
        cut -f 6 {input.SSU} {input.LSU} > {output.summary}
        mv -t {output.wd} {input}
        """
