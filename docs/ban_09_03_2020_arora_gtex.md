## Notes for analysis proposed for Ban (9/3/2020)

### Differences on pipeline processing on GTEx/TCGA data 

**Aims:** 

1. To recreate Arora et al. clustering by PCA and assess whether clustering is better after normalizing by transcriptome diversity.
2. To assess if there are transcritptome diversity differences across the different processing pipelines outlined in Arora et al.

### Files to use (PCA)

The information of the expression matrices and their location is in the following table (see description here):
`/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/expression_datasets.txt`

Transcrriptome diverisity files can be found here:
`/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/original_sources/`

**Script used for PCA is here:**
`/home/users/paedugar/scripts/FraserLab/transcriptome_diveristy/pipelines/R/`

This script:

- takes all arora expression matrices, calculates transcriptome diverisy (**does not normalize it** -- range is [0,log10(nGenes)] ).
- Normalizes gene expression by simple division
- Performs PCA in both, original gene expression matrices and the normalized one
- Writes PC results to disk

Outputs:

They are located here

`/scratch/users/paedugar/transcriptome_diversity/PCA_arora`

There are two outputs per project, expression estimate type:

1. `{project}_{expressionEstimate}.pca.txt` - tab-separated file, with the values of the first 6 PCs and extra info for each sample. The most important column is **pca_type**, its values are "original" for non-normalized matrix, and "transcriptome_diversity_controlled" for the normalized matrix.

2. `{project}_{expressionEstimate}.Rds` - R object containing the full results of the PCA (prcomp function). It's, the slot "pca" contains the results of the non-normalized matrix, and the slot "pca_controlled_transcriptome_diversity" contains the results of the normalized matrix.

