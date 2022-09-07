checkpoint plasmid_reconstruction:
    input:
        "results/trimmed/{isolate}.1.phix.fastq",
        "results/trimmed/{isolate}.2.phix.fastq",
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
        "results/plasmid_reconstruction/{isolate}/run_spades.yaml",
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
        temp("results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe.bam"),
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
    threads: 0.5 * config["threads"]
    wrapper:
        "v1.2.0/bio/samtools/view"


rule sort_primary_alignment:
    input:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.bam",
    output:
        "results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam",
    log:
        "logs/plasmid_reconstruction/{isolate}_sort_primary_aligment.log",
    threads: 0.5 * config["threads"]
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


def guess_longest_kmer(wildcards):
    # Fetch the automatically generated spades configuration file
    spades_conf = checkpoints.plasmid_reconstruction.get(**wildcards).output[1]
    with spades_conf.open() as spades_file:
        # Load as a YAML each of the steps of the analysis
        spades_yaml = yaml.safe_load(spades_file)

    # Extract the steps named after the kmer lengths: K21, K33, K55 etc.
    kmers = [step["STAGE"] for step in spades_yaml if step["STAGE"][0] == "K"]
    # Remove the first character and convert the string to integer
    lengths = [int(x[1:]) for x in kmers]
    # Sort descending order
    lengths.sort(reverse=True)
    # Get the longest
    return lengths[0]


checkpoint plasmid_extraction:
    input:
        graph="results/plasmid_reconstruction/{isolate}/assembly_graph.fastg",
        bam="results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam",
        index="results/plasmid_reconstruction/{isolate}/{isolate}_reads_pe_primary.sorted.bam.bai",
    output:
        "results/plasmid_reconstruction/{isolate}/assembly_graph.cycs.fasta",
    log:
        "logs/plasmid_reconstruction/{isolate}_plasmid_extraction.log",
    params:
        max_kmer_length=guess_longest_kmer,
    conda:
        "../envs/recycler.yaml"
    threads: config["threads"]
    shell:
        # The || solution is from https://stackoverflow.com/a/59055330
        """
        recycle.py -g {input.graph} -b {input.bam} \
            -k {params.max_kmer_length} &> {log} || \
        touch {output} && echo "Empty file produced because of Recycler error" >> {log}
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
            temp("results/trimmed/{isolate}.1.phix.noplasmid.fastq"),
            temp("results/trimmed/{isolate}.2.phix.noplasmid.fastq"),
        ],
    log:
        "logs/plasmid_reconstruction/{isolate}_remove_plasmid_from_reads.log",
    params:
        extra=lambda w, input: "ref={} k=31 hdist=1".format(input.adapters),
    threads: config["threads"]
    wrapper:
        "v1.1.0/bio/bbtools/bbduk"
