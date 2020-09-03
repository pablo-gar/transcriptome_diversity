## Important locations

### Scripts

All scripts are stored in Sherlock at:

`/home/users/paedugar/scripts/FraserLab/transcriptome_diveristy/`

The structure of the folder is simple:

- `./R/` contains R functions used across many other scripts.
- `./fastqtl/` local installation of [fastql](http://fastqtl.sourceforge.net/) used for eQTL searches. 
- `./pipelines/` contains all the pipelines and necessary scripts for data downloading, processing, analyses, and making figures.
   - `./pipelines/R/` all R scripts necessary for downloading and processing data.
   - `./pipelines/figures/` scripts and data to create figures, also the figures themselves. This will likely a stand-alone location to recreate all the figures in the paper.
   - `./pipelines/analysis/` quick and dirty notebook-style scripts to analyze data.
   

### Data

All data is currently stored in Sherlock. However all pipelines are designed to easily modify this for reproducible results.

The most important locations are defined in `config.json` located in the home scripts folder (see above). This file is used in snakemake pipelines as a reference to identify where to store data and results.

The important slots in `config.json` are (this is not comprehensive):

- `projectDir` this is the root folder, **all other folders are relative to it.** Currently set as `/scratch/users/paedugar/transcriptome_diversity`
- `vcfFileGTEXCommonV8` contains the full path to the GTEx genotype VCF file
- `expression_datasets_dir` relative path to data downloaded from orignal sources.
- `expression_matrices_dir*` relative path to folder containing standardized expression matrices
- `transcriptome_diversity_dir` relative path to folder containing transcriptome diversity values
- `auxiliary_Files` list of auxiliary items:
    - `dir` root dir.
    - `expression_matrix_files` relative path to a table describing all expression matrices (TPM, UP-quant, RPKM) as downloaded from **original sources.**
    - `expression_matrix_files_gtex_bed` relative path to a table describing expression matrices from **GTEx in bed format** (used for eQTL searches).
    - `expression_matrix_files_gtex_tpm` relative path to a table describing expression matrices from **GTEx in TPM format.**
    - `expression_matrix_files_for_expression_associations`  relative path to a table describing expression matrices with **raw counts**.

### Sample tables

From the `config.json` file the following two tables are of outmost importance

#### Expression matrices as described in original sources

Location 

- General:

`[projectDir]/[auxiliary_Files{dir}]/[auxiliary_Files{expression_matrix_files}]`

- Current:

`/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/expression_datasets.txt`

Columns (tab-separated):

1. id - id of expression matrix
2. name - currently same as ID
3. path - location of matrix relative to `[projectDir]` in `config.json`
4. group - original source (e.g. "arora")
5. subgroup - project (e.g. "gtex") 
6. method - processing pipeline (e.g. "mskccBatch)
7. count_type - expression estimate (e.g. "TPM")
8. tissue - tissue/body site/cancer (e.g. "Whole_Blood")

#### Expression matrices with raw counts

Location 

- General:

`[projectDir]/[auxiliary_Files{dir}]/[auxiliary_Files{expression_matrix_files_for_expression_associations}]`

- Current:

`/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/expression_datasets_raw_counts.txt`

Columns (tab-separated):

1. id - id of expression matrix
2. name - currently same as ID
3. path - location of matrix relative to `[projectDir]` in `config.json`
4. group - original source 
5. subgroup - project 
6. method - processing pipeline
7. count_type - expression estimate
8. tissue - tissue/body site/cancer
9. gtf - location of annotation file used for TPM calculations relative to `[projectDir]` in `config.json`

