rule remove_adaptater_filter_length:
    input:
        get_fastqs,
    output:
        "results/{isolate}.fasta",
    shell:
        "head -n 1 {input} > {output}"
