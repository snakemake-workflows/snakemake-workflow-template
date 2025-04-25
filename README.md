# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

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

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Workflow overview

This workflow is a best-practice workflow for `<detailed description>`.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Download genome reference from NCBI
2. Simulate short read sequencing data on the fly (`dwgsim`)
3. Check quality of input read data (`FastQC`)
4. Collect statistics from tool output (`MultiQC`)

## Running the workflow

### Input data

This template workflow creates artificial sequencing data in `*.fastq.gz` format. It does not contain actual input data. The simulated input files are nevertheless created based on a mandatory table linked in the `config.yml` file (default: `.test/samples.tsv`). The sample sheet has the following layout:

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
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

To run the workflow with **apptainer** / **singularity**, add a link to a container registry in the `Snakefile`, for example:
`container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry. Run the workflow with:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

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
| random_reads       | num  | frequency of random read sequences    | 0.01                           |

## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO

- Replace `<owner>` and `<repo>` everywhere in the template with the correct user name/organization, and the repository name. The workflow will be automaticallky added to the [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html) once it is publicly available on Github.
- Replace `<name>` with the workflow name (can be the same as `<repo>`).
- Replace `<description>` with a description of what the workflow does.
- Update the [workflow overview](#running-the-workflow), and [running instructions](#running-the-workflow) including parameters, deployment, authors and references
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment options](#execution)
- Do not forget to also adjust the configuration-specific `config/README.md` file
