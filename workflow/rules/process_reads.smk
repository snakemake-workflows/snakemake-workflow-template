# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


# fetch genome sequence from NCBI
# -----------------------------------------------------
rule get_genome:
    output:
        fasta="results/get_genome/genome.fna",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Downloading genome sequence."""
    params:
        ncbi_ftp=config["get_genome"]["ncbi_ftp"],
    log:
        "results/get_genome/genome.log",
    shell:
        "wget -O results/get_genome/genome.fna.gz {params.ncbi_ftp} > {log} 2>&1 && "
        "gunzip results/get_genome/genome.fna.gz >> {log} 2>&1"


# simulate read data using DWGSIM
# -----------------------------------------------------
rule simulate_reads:
    input:
        fasta=rules.get_genome.output.fasta,
    output:
        fastq1="results/simulate_reads/{sample}.bwa.read1.fastq.gz",
        fastq2="results/simulate_reads/{sample}.bwa.read2.fastq.gz",
    conda:
        "../envs/simulate_reads.yml"
    message:
        """--- Simulating read data with DWGSIM."""
    params:
        read_length=config["simulate_reads"]["read_length"],
        read_number=config["simulate_reads"]["read_number"],
        random_reads=config["simulate_reads"]["random_reads"],
    log:
        "results/simulate_reads/{sample}.log",
    script:
        "../scripts/simulate_reads.py"


# make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="results/simulate_reads/{sample}.bwa.{read}.fastq.gz",
    output:
        html="results/fastqc/{sample}.bwa.{read}_fastqc.html",
        zip="results/fastqc/{sample}.bwa.{read}_fastqc.zip",
    params:
        extra="--quiet",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "results/fastqc/{sample}.bwa.{read}.log",
    threads: 1
    wrapper:
        "v6.0.0/bio/fastqc"


# run multiQC on tool output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}.bwa.{read}_fastqc.html",
            sample=samples.index,
            read=["read1", "read2"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra="--verbose --dirs",
    message:
        """--- Generating MultiQC report for seq data."""
    log:
        "results/multiqc/multiqc.log",
    wrapper:
        "v6.0.0/bio/multiqc"
