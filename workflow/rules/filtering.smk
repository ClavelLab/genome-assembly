rule remove_adaptater_filter_length:
    input:
        unpack(get_fastqs),
    output:
        r1="results/trimmed/{isolate}.1.fastq",
        r2="results/trimmed/{isolate}.2.fastq",
        r1_unpaired="results/trimmed/{isolate}.1.unpaired.fastq",
        r2_unpaired="results/trimmed/{isolate}.2.unpaired.fastq",
    log:
        "logs/remove_adaptater/{isolate}.log",
    conda:
        "envs/trimmomatic.yaml"
    params:
        illumina_clip=config["adapters"],
        extra="LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50",
    shell:
        """
        trimmomatic PE -phred33 {input.r1} {input.r2} \
        {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} \
        ILLUMINACLIP:{params.illumina_clip}:2:30:10 \
        {params.extra}
        """
