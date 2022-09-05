rule remove_adapters_filter_length:
    input:
        unpack(get_fastqs),
    output:
        r1=temp("results/trimmed/{isolate}.1.fastq"),
        r2=temp("results/trimmed/{isolate}.2.fastq"),
        r1_unpaired=temp("results/trimmed/{isolate}.1.unpaired.fastq"),
        r2_unpaired=temp("results/trimmed/{isolate}.2.unpaired.fastq"),
    log:
        "logs/remove_adapters/{isolate}.log",
    conda:
        "../envs/trimmomatic.yaml"
    threads: config["threads"]
    params:
        illumina_clip=config["adapters"],
        extra="LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50",
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
        {input.r1} {input.r2} \
        {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} \
        ILLUMINACLIP:{params.illumina_clip}:2:30:10 \
        {params.extra} &> {log}
        """


rule remove_phix:
    input:
        sample=[
            "results/trimmed/{isolate}.1.fastq",
            "results/trimmed/{isolate}.2.fastq",
        ],
        adapters=config["phix"],
    output:
        trimmed=[
            temp("results/trimmed/{isolate}.1.phix.fastq"),
            temp("results/trimmed/{isolate}.2.phix.fastq"),
        ],
    log:
        "logs/remove_phix/{isolate}.log",
    params:
        extra=lambda w, input: "ref={} k=31 hdist=1".format(input.adapters),
    threads: config["threads"]
    wrapper:
        "v1.1.0/bio/bbtools/bbduk"
