# ictv-taxonomy-challenge-nr

This repository contains the code necessary to reproduce the `ictv-taxonomy-challenge-nr` submission to the 
[ICTV Computational Virus Taxonomy Challenge](https://ictv-vbeg.github.io/ICTV-TaxonomyChallenge/).

## Methods

This submission utilizes a pipeline in which [`mmseqs2`](https://github.com/soedinglab/MMseqs2) is employed to identify ORFs within the input sequences and classify them into ICTV taxa by aligning them against a [reference database](https://doi.org/10.5281/zenodo.14889287) of viral proteins. The protein sequences in this database were sourced from NCBI's NR database on 2025-02-08 and mapped to a custom taxonomy generated using [`taxonkit`](https://github.com/shenwei356/taxonkit).

The code in this repository enables the reproducibility of genome classification in the ICTV Computational Virus Taxonomy Challenge. However, it does not include the scripts required to build the custom protein database. For that, please refer to the [`apcamargo/ictv-mmseqs2-protein-database`](https://github.com/apcamargo/ictv-mmseqs2-protein-database) repository. Note that generating a database from scratch will result in a dataset different from the one used in this submission, as the NR database is continuously updated.

## Reproduction

This repository is structured as a [Pixi](https://pixi.sh/) project, which facilitates the recreation of the environment with the correct dependencies and enables execution of the code. Once you've cloned the repository, you can run the pipeline by executing the following command:

```sh
pixi run pipeline
```

## Authors
- [@apcamargo](https://github.com/apcamargo)
- [@UriNeri](https://github.com/UriNeri)
- [@LanderDC](https://github.com/LanderDC)
