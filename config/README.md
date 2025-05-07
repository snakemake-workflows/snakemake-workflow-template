## Workflow overview

This workflow is a best-practice workflow for `<detailed description>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Download genome reference from NCBI
2. Validate downloaded genome (`python` script)
3. Simulate short read sequencing data on the fly (`dwgsim`)
4. Check quality of input read data (`FastQC`)
5. Collect statistics from tool output (`MultiQC`)

## Running the workflow

### Input data

This template workflow creates artificial sequencing data in `*.fastq.gz` format.
It does not contain actual input data.
The simulated input files are nevertheless created based on a mandatory table linked in the `config.yml` file (default: `.test/samples.tsv`).
The sample sheet has the following layout:

| sample  | condition | replicate | read1                      | read2                      |
| ------- | --------- | --------- | -------------------------- | -------------------------- |
| sample1 | wild_type | 1         | sample1.bwa.read1.fastq.gz | sample1.bwa.read2.fastq.gz |
| sample2 | wild_type | 2         | sample2.bwa.read1.fastq.gz | sample2.bwa.read2.fastq.gz |

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter          | type | details                               | default                        |
| ------------------ | ---- | ------------------------------------- | ------------------------------ |
| **samplesheet**    |      |                                       |                                |
| path               | str  | path to samplesheet, mandatory        | "config/samples.tsv"           |
| **get_genome**     |      |                                       |                                |
| ncbi_ftp           | str  | link to a genome on NCBI's FTP server | link to _S. cerevisiae_ genome |
| **simulate_reads** |      |                                       |                                |
| read_length        | num  | length of target reads in bp          | 100                            |
| read_number        | num  | number of total reads to be simulated | 10000                          |
