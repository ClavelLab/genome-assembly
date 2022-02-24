rule plasmid_reconstruction:
    input:
        "results/trimmed/{isolate}.1.phix.fastq",
        "results/trimmed/{isolate}.2.phix.fastq",
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    log:
        "logs/plasmid_reconstruction/{isolate}_plasmid_reconstruction.log",
    conda:
        "../envs/spades.yaml"
    threads: config["threads"]
    shell:
        """
        plasmidspades.py -t {threads} --careful \
            -1 {input[0]} -2 {input[1]} \
            -o results/plasmid_reconstruction/{wildcards.isolate} &> {log}
        """


rule extract_plasmid_assembly_graph:
    input:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
    log:
        "logs/plasmid_reconstruction/{isolate}_assembly_extraction.log",
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        make_fasta_from_fastg.py -g {input} -o {output} &> {log}
        """


rule index_plasmid_assembly_graph:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
    output:
        idx=multiext(
            "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/plasmid_reconstruction/{isolate}_index_assembly_graph.log",
    wrapper:
        "v1.2.0/bio/bwa/index"


rule align_plasmid_assembly_graph:
    input:
        reads=[
            "results/trimmed/{isolate}.1.phix.fastq",
            "results/trimmed/{isolate}.2.phix.fastq",
        ],
        idx=multiext(
            "results/plasmid_reconstruction/{isolate}/{isolate}_assembly_graph.nodes.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe.bam",
    log:
        "logs/plasmid_reconstruction/{isolate}_align_assembly_graph.log",
    params:
        # No sorting. But defaults include header and compressed BAM
        sorting="none",
    threads: config["threads"]
    wrapper:
        "v1.2.0/bio/bwa/mem"


rule select_primary_alignment:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe.bam",
    output:
        bam=pipe(
            "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.bam"
        ),
    log:
        "logs/plasmid_reconstruction/{isolate}_select_primary_aligment.log",
    params:
        extra="-bF 0x0800",  # remove supplementary alignments
    threads: config["threads"]
    wrapper:
        "v1.2.0/bio/samtools/view"


rule sort_primary_alignment:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.bam",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam",
    log:
        "logs/plasmid_reconstruction/{isolate}_sort_primary_aligment.log",
    threads: config["threads"]
    wrapper:
        "v1.2.0/bio/samtools/sort"


rule index_primary_alignment:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam.bai",
    log:
        "logs/plasmid_reconstruction/{isolate}_index_primary_aligment.log",
    threads: config["threads"]
    wrapper:
        "v1.2.0/bio/samtools/index"


rule plasmid_extraction:
    input:
        graph="results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
        bam="results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam",
        index="results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam.bai",
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    log:
        "logs/plasmid_reconstruction/{isolate}_plasmid_extraction.log",
    params:
        max_kmer_length=55,
    conda:
        "../envs/recycler.yaml"
    threads: config["threads"]
    shell:
        """
        recycle.py -g {input.graph} -b {input.bam} \
            -k {params.max_kmer_length} &> {log}
        """


rule remove_plasmid_from_reads:
    input:
        sample=[
            "results/trimmed/{isolate}.1.phix.fastq",
            "results/trimmed/{isolate}.2.phix.fastq",
        ],
        adapters="results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    output:
        trimmed=[
            "results/trimmed/{isolate}.1.phix.noplasmid.fastq",
            "results/trimmed/{isolate}.2.phix.noplasmid.fastq",
        ],
    log:
        "logs/plasmid_reconstruction/{isolate}_remove_plasmid_from_reads.log",
    params:
        extra=lambda w, input: "ref={} k=31 hdist=1".format(input.adapters),
    threads: config["threads"]
    wrapper:
        "v1.1.0/bio/bbtools/bbduk"
