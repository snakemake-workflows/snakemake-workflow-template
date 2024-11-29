# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MPUSP/snakemake-workflow-template/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/MPUSP/snakemake-workflow-template/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1D355C.svg?labelColor=000000)](https://sylabs.io/docs/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog)

A Snakemake workflow for `<description>`

- [Snakemake workflow: `<name>`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Workflow overview](#workflow-overview)
  - [Running the workflow](#running-the-workflow)
    - [Input data](#input-data)
    - [Execution](#execution)
    - [Parameters](#parameters)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Workflow overview

This workflow is a best-practice workflow for `<detailed description>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Parse sample sheet containing sample meta data (`python`)
2. Simulate short read sequencing data on the fly (`dwgsim`)
3. Check quality of input read data (`FastQC`)
4. Trim adapters from input data (`cutadapt`)
5. Collect statistics from tool output (`MultiQC`)

## Running the workflow

### Input data

This template workflow contains artifical sequencing data in `*.fastq.gz` format.
The test data is located in `.test/data`. Input files are supplied with a mandatory table, whose location is indicated in the `config.yml` file (default: `.test/samples.tsv`). The sample sheet has the following layout:

| sample  | condition | replicate | read1                      | read2                      |
| ------- | --------- | --------- | -------------------------- | -------------------------- |
| sample1 | wild_type | 1         | sample1.bwa.read1.fastq.gz | sample1.bwa.read2.fastq.gz |
| sample2 | wild_type | 2         | sample2.bwa.read1.fastq.gz | sample2.bwa.read2.fastq.gz |


### Execution

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-workflow-name
```

Adjust options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the complete workflow with test files using **conda**, execute the following command. The definition of the number of compute cores is mandatory.

```bash
snakemake --cores 10 --sdm conda --directory .test
```

To run the workflow with **singularity** / **apptainer**, use:

```bash
snakemake --cores 10 --sdm conda apptainer --directory .test
```

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter          | type | details                                 | default                                       |
| ------------------ | ---- | --------------------------------------- | --------------------------------------------- |
| **samplesheet**    |      |                                         |                                               |
| path               | str  | path to samplesheet, mandatory          | "config/samples.tsv"                          |
| **get_genome**     |      |                                         |                                               |
| database           | str  | one of `manual`, `ncbi`                 | `ncbi`                                        |
| assembly           | str  | RefSeq ID                               | `GCF_000006785.2`                             |
| fasta              | str  | optional path to fasta file             | Null                                          |
| gff                | str  | optional path to gff file               | Null                                          |
| gff_source_type    | str  | list of name/value pairs for GFF source | see config file                               |
| **simulate_reads** |      |                                         |                                               |
| read_length        | num  |                                         | 100                                           |
| read_number        |      |                                         | 100000                                        |
| random_freq        |      |                                         | 0.01                                          |
| **cutadapt**       |      |                                         |                                               |
| threep_adapter     | str  | sequence of the 3' adapter              | `-a ATCGTAGATCGG`                             |
| fivep_adapter      | str  | sequence of the 5' adapter              | `-A GATGGCGATAGG`                             |
| default            | str  | additional options passed to `cutadapt` | [`-q 10 `, `-m 25 `, `-M 100`, `--overlap=5`] |
| **multiqc**        |      |                                         |                                               |
| config             | str  | path to multiqc config                  | `config/multiqc_config.yml`                   |

## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. *Sustainable data analysis with Snakemake*. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* Update the workflow description, parameters, running options, authors and references in the `README.md`
* Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's capability
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.