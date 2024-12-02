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
Before running the entire workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the complete workflow with test files using **conda**, execute the following command. The definition of the number of compute cores is mandatory.

```bash
snakemake --cores 3 --sdm conda --directory .test
```

To run the workflow with **singularity** / **apptainer**, add a link to a container registry in the `Snakefile`, for example:
`container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry. Run the workflow with:

```bash
snakemake --cores 3 --sdm conda apptainer --directory .test
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
| read_length        | num  | length of target reads in bp            | 100                                           |
| read_number        | num  | number of total reads to be simulated   | 100000                                        |
| random_freq        | num  | frequency of random read sequences      | 0.01                                          |
| **cutadapt**       |      |                                         |                                               |
| threep_adapter     | str  | sequence of the 3' adapter              | `-a ATCGTAGATCGG`                             |
| fivep_adapter      | str  | sequence of the 5' adapter              | `-A GATGGCGATAGG`                             |
| default            | str  | additional options passed to `cutadapt` | [`-q 10 `, `-m 25 `, `-M 100`, `--overlap=5`] |
| **multiqc**        |      |                                         |                                               |
| config             | str  | path to multiQC config                  | `config/multiqc_config.yml`                   |

## TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* Update the workflow parameters and running options
