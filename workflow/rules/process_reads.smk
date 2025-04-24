# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


# fetch genome from NCBI
# -----------------------------------------------------
rule get_genome:
    output:
        fasta="results/get_genome/genome.fna",
        gff="results/get_genome/genome.gff",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Downloading genome FASTA and GFF files."""
    params:
        ncbi_ftp=config["get_genome"]["ncbi_ftp"],
    log:
        path="results/get_genome/genome.log",
    shell:
        "for file_type in fna gff; do "
        "  wget -O results/get_genome/genome.${{file_type}}.gz "
        "  {params.ncbi_ftp}.${{file_type}}.gz && "
        "  gunzip results/get_genome/genome.${{file_type}}.gz; "
        "done > {log.path} 2>&1"


# simulate reads using DWGSIM
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
        """--- Simulating reads with DWGSIM."""
    params:
        read_length=config["simulate_reads"]["read_length"],
        read_number=config["simulate_reads"]["read_number"],
        random_freq=config["simulate_reads"]["random_freq"],
    log:
        path="results/simulate_reads/{sample}.log",
    shell:
        "prefix=`echo {output.fastq1} | cut -f 1 -d .`;"
        "dwgsim -1 {params.read_length} -2 {params.read_length} "
        "-N {params.read_number} -o 1 -y {params.random_freq} {input.fasta} ${{prefix}} &> {log.path}"


# make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="results/simulate_reads/{sample}.bwa.{read}.fastq.gz",
    output:
        report="results/fastqc/{sample}.bwa.{read}_fastqc.html",
    conda:
        "../envs/fastqc.yml"
    message:
        """--- Checking fastq files with FastQC."""
    log:
        path="results/fastqc/{sample}.bwa.{read}.log",
    shell:
        "prefix=`echo {output.report} | cut -f 1-2 -d /`;"
        "fastqc --nogroup --extract --quiet --threads {threads} -o ${{prefix}} {input.fastq} > {log}"


# trim adapters from reads
# -----------------------------------------------------
rule cutadapt:
    input:
        fastq1="results/simulate_reads/{sample}.bwa.read1.fastq.gz",
        fastq2="results/simulate_reads/{sample}.bwa.read2.fastq.gz",
    output:
        fastq1="results/cutadapt/{sample}.bwa.read1.fastq.gz",
        fastq2="results/cutadapt/{sample}.bwa.read2.fastq.gz",
    conda:
        "../envs/cutadapt.yml"
    message:
        """--- Trim adapters from reads."""
    params:
        threep_adapter=config["cutadapt"]["threep_adapter"],
        fivep_adapter=config["cutadapt"]["fivep_adapter"],
        default=config["cutadapt"]["default"],
    log:
        stdout="results/cutadapt/{sample}.bwa.log",
        stderr="results/cutadapt/{sample}.bwa.stderr",
    threads: int(workflow.cores * 0.25)
    shell:
        "cutadapt {params.threep_adapter} {params.fivep_adapter} --cores {threads} "
        "-o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq1} > {log.stdout} 2> {log.stderr}"


# run multiQC on tool output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/cutadapt/{sample}.bwa.{read}.fastq.gz",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        expand(
            "results/fastqc/{sample}.bwa.{read}_fastqc.html",
            sample=samples.index,
            read=["read1", "read2"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    message:
        """--- Generating MultiQC report for seq data."""
    params:
        config=config["multiqc"]["config"],
    log:
        path="results/multiqc/log/multiqc.log",
    shell:
        "outdir=`echo {output.report} | cut -f 1-2 -d /`; "
        "multiqc -c {params.config} --force --verbose --dirs --outdir ${{outdir}} results &> {log.path}"
