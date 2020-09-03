# Notes for analysis proposed for Ban (9/3/2020)


## Differences on pipeline processing on GTEx/TCGA data 

###Aims

1. To recreate Arora et al. clustering by PCA and assess whether clustering is better after normalizing by transcriptome diversity.
2. To assess if there are transcritptome diversity differences across the different processing pipelines outlined in Arora et al.

### Files to use

The information of the expression matrices and their location is in the following table (see description [here](https://github.com/pablo-gar/transcriptome_diversity/blob/master/docs/important_locations.md)):
`/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/expression_datasets.txt`

Transcrriptome diversity files can be found here:
`/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/original_sources/`

**Script used for PCA is here:**
`/home/users/paedugar/scripts/FraserLab/transcriptome_diveristy/pipelines/R/`

This script:

- takes all arora expression matrices, calculates transcriptome diversify (**does not normalize it** -- range is [0,log10(nGenes)] ).
- Normalizes gene expression by simple division
- Performs PCA in both, original gene expression matrices and the normalized one
- Writes PC results to disk

Outputs:

They are located here

`/scratch/users/paedugar/transcriptome_diversity/PCA_arora`

There are two outputs for each project/expression-estimate combination:

1. `{project}_{expressionEstimate}.pca.txt` - tab-separated file, with the values of the first 6 PCs and extra info for each sample. The most important column is **pca_type**, its values are "original" for non-normalized matrix, and "transcriptome_diversity_controlled" for the normalized matrix.

2. `{project}_{expressionEstimate}.Rds` - R object containing the full results of the PCA (prcomp function), as well as the original matrices. It's a list, the slot "pca" contains the results of the non-normalized matrix, and the slot "pca_controlled_transcriptome_diversity" contains the results of the normalized matrix.

## Associations between GTEx metadata and transcriptome diversity

### Aims:

1. To assess if and which metadata information (sex, age, RIN ischemic time, etc.) are associated with the transcriptome diversity of each sample
2. To assess if there are significant differences on transcriptome diversity across tissues

### Files to use

Transcriptome diversity files are here
`/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm`

GTEx metadata are here:
`/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/GTEx_Phenotypes_v8`

The two most important files in that folder are:
`phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz` Contains the metadata per person
`phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz` Contains the metadata per sample

