# Snakemake workflow: `genome-assembly-clavel`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/cpauvert/genome-assembly-clavel/workflows/Tests/badge.svg?branch=main)](https://github.com/cpauvert/genome-assembly-clavel/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `assembling bacterial genomes according to the standard operating procedure in the Clavel Lab`


## Usage

First, the conda environment can be installed using the following command

    snakemake --use-conda --conda-create-envs-only

Then, the workflow can be used with the default configuration file and 9 threads here:

    snakemake -c 9 --use-conda

Or using a custom one:

    snakemake --configfile config/er2.yaml -c 9 --use-conda

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cpauvert%2Fgenome-assembly-clavel).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository its DOI (see above).

