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
        ncbi_ftp=lookup(within=config, dpath="get_genome/ncbi_ftp"),
    log:
        "results/get_genome/genome.log",
    shell:
        "wget -O results/get_genome/genome.fna.gz {params.ncbi_ftp} > {log} 2>&1 && "
        "gunzip results/get_genome/genome.fna.gz >> {log} 2>&1"


# validate genome sequence file
# -----------------------------------------------------
rule validate_genome:
    input:
        fasta=rules.get_genome.output.fasta,
    output:
        fasta="results/validate_genome/genome.fna",
    conda:
        "../envs/validate_genome.yml"
    message:
        """--- Validating genome sequence file."""
    log:
        "results/validate_genome/genome.log",
    script:
        "../scripts/validate_fasta.py"


# simulate read data using DWGSIM
# -----------------------------------------------------
rule simulate_reads:
    input:
        fasta=rules.validate_genome.output.fasta,
    output:
        multiext(
            "results/simulate_reads/{sample}",
            read1=".bwa.read1.fastq.gz",
            read2=".bwa.read2.fastq.gz",
        ),
    conda:
        "../envs/simulate_reads.yml"
    message:
        """--- Simulating read data with DWGSIM."""
    params:
        output_type=1,
        read_length=lookup(within=config, dpath="simulate_reads/read_length"),
        read_number=lookup(within=config, dpath="simulate_reads/read_number"),
    log:
        "results/simulate_reads/{sample}.log",
    shell:
        "output_prefix=`echo {output.read1} | cut -f 1 -d .`;"
        "dwgsim "
        " -1 {params.read_length}"
        " -2 {params.read_length}"
        " -N {params.read_number}"
        " -o {params.output_type}"
        " {input.fasta}"
        " ${{output_prefix}}"
        " > {log} 2>&1"


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
            "results/fastqc/{sample}.bwa.{read}_fastqc.{ext}",
            sample=samples.index,
            read=["read1", "read2"],
            ext=["html", "zip"],
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
