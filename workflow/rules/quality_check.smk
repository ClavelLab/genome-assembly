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
    priority: 10
    params:
        min_contig_length=config["min_contig_length"]
    shell:
        """
        seqkit seq --remove-gaps --min-len {params.min_contig_length} --threads {threads} {input} 1> {output} 2> {log}
        """


rule rewrite_genome_headers:
    input:
        "results/quality_check/{isolate}/trimmed/{isolate}.trimmed.fa",
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
        "results/quality_check/{isolate}/bakta/{isolate}.tsv",
    log:
        "logs/quality_check/{isolate}_bakta.log",
    params:
        db=config["bakta_db"],
        outdir="results/quality_check/{isolate}/bakta",
    conda:
        "../envs/bakta.yaml"
    threads: config["threads"]
    shell:
        """
        bakta --db {params.db}/db/ \
        --prefix {wildcards.isolate} \
        --locus-tag {wildcards.isolate} \
        --compliant \
        --output {params.outdir} \
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
        "results/quality_check/{isolate}/checkm/{isolate}_checkm.tsv",
        seq="SSU-{isolate}.extraction.fasta",
    output:
        "results/quality_check/{isolate}/{isolate}.SSU.fa",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        seqkit replace -p 'NODE' -r '{wildcards.isolate}_SSU_16S_rRNA NODE' --line-width 0 --threads {threads} {input.seq} > {output}
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
        "results/quality_check/{isolate}/checkm/{isolate}_checkm.tsv",
        seq="LSU-{isolate}.extraction.fasta",
    output:
        "results/quality_check/{isolate}/{isolate}.LSU.fa",
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    shell:
        """
        seqkit replace -p 'NODE' -r '{wildcards.isolate}_LSU_23S_rRNA NODE' --line-width 0 --threads {threads} {input.seq} > {output}
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


use rule compute_lengths_plasmids as compute_lengths_SSU with:
    input:
        "results/quality_check/{isolate}/{isolate}.SSU.fa",
    output:
        "results/quality_check/{isolate}/metaxa/{isolate}.SSU.tsv",
    log:
        "logs/quality_check/{isolate}_compute_lengths_SSU.log",


use rule compute_lengths_plasmids as compute_lengths_LSU with:
    input:
        "results/quality_check/{isolate}/{isolate}.LSU.fa",
    output:
        "results/quality_check/{isolate}/metaxa/{isolate}.LSU.tsv",
    log:
        "logs/quality_check/{isolate}_compute_lengths_LSU.log",


rule aggregate_SSU_LSU_results:
    input:
        SSU="results/quality_check/{isolate}/metaxa/{isolate}.SSU.tsv",
        LSU="results/quality_check/{isolate}/metaxa/{isolate}.LSU.tsv",
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
    priority: 30
    shell:
        """
        seqkit stats -T --threads {threads} {input} 1> {output} 2> {log}
        """


use rule basepairs_metrics as trimmed_basepairs_metrics with:
    input:
        "results/trimmed/{isolate}.1.fastq",
        "results/trimmed/{isolate}.2.fastq",
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.trimmed.tsv",
    log:
        "logs/quality_check/{isolate}_trimmed_basepairs_metrics.log",


use rule basepairs_metrics as assembly_basepairs_metrics with:
    input:
        "results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        "results/quality_check/{isolate}/metrics/{isolate}.assembly.tsv",
    log:
        "logs/quality_check/{isolate}_assembly_basepairs_metrics.log",


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
    priority: 50
    script:
        "../scripts/write_coverage_and_metrics.py"


rule quast_for_assembly_quality:
    input:
        raw_assembly="results/assembly/{isolate}/contigs.fasta",
        final_assembly="results/quality_check/{isolate}/{isolate}.genome.fa",
    output:
        quast_dir=directory("results/quality_check/{isolate}/quast"),
        report="results/quality_check/{isolate}/quast/transposed_report.tsv",
    log:
        "logs/quality_check/{isolate}_quast.log",
    conda:
        "../envs/quast.yaml"
    threads: config["threads"]
    shell:
        """
        quast --threads {threads} --labels "{wildcards.isolate}.raw,{wildcards.isolate}.final" \
        --no-icarus --min-contig 0 --output-dir {output.quast_dir} {input.raw_assembly} {input.final_assembly} > {log}
        """


rule checksum_raw_fastq:
    input:
        unpack(get_fastqs),
    output:
        "results/quality_check/{isolate}/checksums/{isolate}_raw_fastq.md5",
    log:
        "logs/quality_check/{isolate}_checksum_fastq.log",
    shell:
        "md5sum {input} 1> {output} 2> {log}"


rule checksum_genome:
    input:
        "results/genome/{isolate}.genome.fa.gz",
    output:
        "results/quality_check/{isolate}/checksums/{isolate}_genome.md5",
    log:
        "logs/quality_check/{isolate}_checksum_genome.log",
    shell:
        "md5sum {input} 1> {output} 2> {log}"


rule write_summary_table:
    input:
        samples=config["samples"],
        metrics="results/quality_check/{isolate}/metrics/{isolate}.metrics.csv",
        fastq_md5="results/quality_check/{isolate}/checksums/{isolate}_raw_fastq.md5",
        checkm="results/quality_check/{isolate}/checkm/{isolate}_checkm.tsv",
        ssu_lsu="results/quality_check/{isolate}/metaxa/{isolate}.SSU-LSU.csv",
        trnas_5s="results/quality_check/{isolate}/bakta/{isolate}.tRNAs-5S.csv",
        quast="results/quality_check/{isolate}/quast/transposed_report.tsv",
        genome_md5="results/quality_check/{isolate}/checksums/{isolate}_genome.md5",
        plasmids="results/plasmid_reconstruction/{isolate}/plasmids_lengths.tsv",
    output:
        "results/summary/{isolate}.csv",
    log:
        "logs/quality_check/{isolate}_summary_table.log",
    params:
        adapters=config["adapters"],
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/write_summary_table.py"


rule aggregate_summaries:
    input:
        expand(
            "results/summary/{isolate}.csv",
            isolate=samples["isolate"],
        ),
    output:
        "results/"
        + os.path.basename(config["samples"]).removesuffix(".tsv")
        + "-summary.csv",
    log:
        "logs/quality_check/aggregate_summaries.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/aggregate_summaries.py"
