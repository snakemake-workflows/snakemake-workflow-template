# A best practice template for Snakemake workflows

To generate a new structure using this template, run

    copier gh:snakemake-workflows/snakemake-workflow-template <path>

with path pointing to the desired output directory.

You will get the following structure:

    ├── copier.yml
    ├── README.md
    ├── .template
    │   ├── config
    │   │   └── config.yaml.tmpl
    │   └── workflow
    │       └── Snakefile.tmpl
    └── workflow
        ├── rules
        │   └── common.smk
        └── Snakefile

The `workflow` folder contains the structure to put the actual workflow in.
Make sure to follow the [best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility).

The `.template` folder contains an automatically generated template for deploying the workflow using Snakemake's [module system](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules).
Here, you should edit the `config.yaml.tmpl` to contain the config settings you would like to initialize in the deployment of the user.
It is possible to add further files needed for the configuration to the `.template` folder, e.g. sample sheets.