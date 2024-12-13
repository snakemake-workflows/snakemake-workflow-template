# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


# module to fetch genome from NCBI or Ensemble
# -----------------------------------------------------
rule get_genome:
    output:
        path=directory("results/get_genome"),
        fasta="results/get_genome/genome.fasta",
        gff="results/get_genome/genome.gff",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Parsing genome GFF and FASTA files."""
    params:
        database=config["get_genome"]["database"],
        assembly=config["get_genome"]["assembly"],
        fasta=config["get_genome"]["fasta"],
        gff=config["get_genome"]["gff"],
    log:
        path="results/get_genome/log/get_genome.log",
    script:
        "../scripts/get_genome.py"


# module simulate reads
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


# module to make QC report
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


# module to trim adapters from reads
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


# module to run multiQC on input + processed files
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
