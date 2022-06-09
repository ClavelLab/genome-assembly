rule remove_small_contigs:
    input:
        "results/assembly/{isolate}/contigs.fasta",
    output:
        "results/quality_check/{isolate}/trimmed/{isolate}.trimmed.fa",
    log:
        "logs/quality_check/{isolate}_remove_small_contigs.log",
    conda:
        "../envs/seqkit.yaml"
    threads: config["threads"]
    shell:
        """
        seqkit seq --remove-gaps --min-len 1000 --threads {threads} {input} 1> {output} 2> {log}
        """


rule download_mdmcleaner_db:
    input:
        HTTP.remote(
            "zenodo.org/record/5698995/files/MDMcleanerDB.tar.bz2", keep_local=False
        ),
    output:
        directory(config["mdmcleaner_db"]),
    log:
        "logs/quality_check/download_mdmcleaner_db.log",
    conda:
        "../envs/mdmcleaner.yaml"
    shell:
        """
        mkdir -p {output}
        # change directory before extracting
        tar -C {output} -xvjf {input} 2>&1 {log}
        """


rule set_mdmcleaner_db:
    output:
        "results/quality_check/mdmcleaner.config",
    log:
        "logs/quality_check/set_mdmcleaner_db.log",
    params:
        db=config["mdmcleaner_db"],
    conda:
        "../envs/mdmcleaner.yaml"
    shell:
        """
        mdmcleaner set_configs --db_basedir {params.db}
        mv mdmcleaner.config {output}
        """


rule mdmcleaner_for_contamination_check:
    input:
        genome="results/quality_check/{isolate}/trimmed/{isolate}.trimmed.fa",
        config="results/quality_check/mdmcleaner.config",
    output:
        mdm_dir=directory("results/quality_check/{isolate}/mdmcleaner/"),
        filtered_contigs="results/quality_check/{isolate}/mdmcleaner/{isolate}.trimmed_filtered_kept_contigs.fasta.gz",
    log:
        "logs/quality_check/{isolate}_mdmcleaner.log",
    params:
        local_dir="mdmcleaner_{isolate}",
    conda:
        "../envs/mdmcleaner.yaml"
    threads: config["threads"]
    shell:
        """
        mdmcleaner clean --config {input.config} \
        --input_fastas {input.genome} \
        --output_folder {params.local_dir} \
        --threads {threads} &> {log}
        mv {params.local_dir} {output.mdm_dir}
        """


rule rewrite_genome_headers:
    input:
        "results/quality_check/{isolate}/mdmcleaner/{isolate}.trimmed_filtered_kept_contigs.fasta.gz",
    output:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        seqkit replace --pattern 'NODE(_[0-9]*_length_[0-9]*)_cov.*' \
        --replacement '{wildcards.isolate}_contig$1' \
        {input} > {output}
        """


rule checkM_for_quality:
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        "results/quality_check/{isolate}/checkm/{isolate}_checkm.tsv",
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
        seqkit replace -p 'NODE' -r '{wildcards.isolate}_SSU_16S_rRNA NODE' --line-width 0 --threads {threads} {input} > {output}
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
        seqkit replace -p 'NODE' -r '{wildcards.isolate}_LSU_23S_rRNA NODE' --line-width 0 --threads {threads} {input} > {output}
        """


rule move_SSU_LSU_results:
    input:
        SSU="SSU-{isolate}.extraction.results",
        LSU="LSU-{isolate}.extraction.results",
    output:
        SSU="results/quality_check/{isolate}/metaxa/SSU-{isolate}.extraction.results",
        LSU="results/quality_check/{isolate}/metaxa/LSU-{isolate}.extraction.results",
    shell:
        """
        mv {input.SSU} {output.SSU}
        mv {input.LSU} {output.LSU}
        """


rule aggregate_SSU_LSU_results:
    input:
        SSU="results/quality_check/{isolate}/{isolate}.SSU.fa",
        LSU="results/quality_check/{isolate}/{isolate}.LSU.fa",
        results_SSU="results/quality_check/{isolate}/metaxa/SSU-{isolate}.extraction.results",
        results_LSU="results/quality_check/{isolate}/metaxa/LSU-{isolate}.extraction.results",
    output:
        "results/quality_check/{isolate}/metaxa/{isolate}.SSU-LSU.csv",
    log:
        "logs/quality_check/{isolate}_aggregate_SSU_LSU.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/check_LSU_SSU.py"


rule basepairs_metrics:
    input:
        unpack(get_fastqs),
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.raw_sequence.tsv",
    log:
        "logs/quality_check/{isolate}_basepairs_metrics.log",
    conda:
        "../envs/seqkit.yaml"
    threads: config["threads"]
    shell:
        """
        seqkit stats -T --threads {threads} {input} 1> {output} 2> {log}
        """


rule trimmed_basepairs_metrics:
    input:
        "results/trimmed/{isolate}.1.fastq",
        "results/trimmed/{isolate}.2.fastq",
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.trimmed.tsv",
    log:
        "logs/quality_check/{isolate}_trimmed_basepairs_metrics.log",
    conda:
        "../envs/seqkit.yaml"
    threads: config["threads"]
    shell:
        """
        seqkit stats -T --threads {threads} {input} 1> {output} 2> {log}
        """


rule assembly_basepairs_metrics:
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.assembly.tsv",
    log:
        "logs/quality_check/{isolate}_assembly_basepairs_metrics.log",
    conda:
        "../envs/seqkit.yaml"
    threads: config["threads"]
    shell:
        """
        seqkit stats -T --threads {threads} {input} 1> {output} 2> {log}
        """


rule write_coverage_and_metrics:
    input:
        raw="results/quality_check/{isolate}/metrics/{isolate}.raw_sequence.tsv",
        trimmed="results/quality_check/{isolate}/metrics/{isolate}.trimmed.tsv",
        assembly="results/quality_check/{isolate}/metrics/{isolate}.assembly.tsv",
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.metrics.csv",
    log:
        "logs/quality_check/{isolate}_write_coverage_and_metrics.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/write_coverage_and_metrics.py"


rule quast_for_assembly_quality:
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        directory("results/quality_check/{isolate}/quast"),
    log:
        "logs/quality_check/{isolate}_quast.log",
    conda:
        "../envs/quast.yaml"
    threads: config["threads"]
    shell:
        """
        quast --threads {threads} --labels "{wildcards.isolate}" \
        --no-icarus --output-dir {output} {input} > {log}
        """


rule checksum_raw_fastq:
    input:
        unpack(get_fastqs),
    output:
        "results/quality_check/{isolate}/checksums/{isolate}_raw_fastq.md5",
    log:
        "logs/quality_check/{isolate}_checksum_fastq.log",
    shell:
        """
        md5sum {input} 1> {output} 2> {log}
        """
