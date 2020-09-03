# To make expression tables first run:
# cd R; Rscript create_sample_table_DNA_samples.R

import pandas as pd
import os

configfile: '../config.json'
localrules: all, get_chromosome_arms, merge_count_reads_chromosome_arms

keele_tissues=['kidney', 'liver', 'lung']
#keele_tissues=['kidney']

#---------------------------------#
# Load tissue sample table
sample_table = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files']))
sample_table = sample_table.set_index('id')

all_datasets = set(sample_table.index.get_level_values('id'))

#---------------------------------#
# Pipeline
# 

rule all:
    input: 
        # Get transcriptome diversity
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', '{dataset}.txt'), dataset=list(all_datasets)),
        ## Get associations with gene expression 
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations', '{dataset}.txt'), dataset=list(all_datasets)),
        ### Do analysis of PCs and transcriptome diversity
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'PC_matrix.txt'), dataset=list(all_datasets)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'PC_stats.txt'), dataset=list(all_datasets)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'cors.txt'), dataset=list(all_datasets)),
        #### Do PCA analysis on multiple processing pipelines, Arora et al
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_arora', '{project}_{count_type}.Rds'), project=['gtex', 'TCGA'], count_type=['TPM', 'RPKM']),
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_arora', '{project}_{count_type}.pca.txt'), project=['gtex', 'TCGA'], count_type=['TPM', 'RPKM'])
        #
        ##Do differential_gene_expression
        #os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_genotype.Rds'), 
        #os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_fulvestrant_wt.Rds'), 
        #os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_tamoxifen_wt.Rds'), 
        #os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_jq1_wt.Rds'),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'sarantopoulou_{platform}', 'edgeR_results_ILB.Rds'),  platform=['truseq', 'v4']),
        #
        ## Do eqtl search on keele data
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_no_covariate', '{keele_tissue}.{chunk}.txt.gz'), keele_tissue=keele_tissues, chunk=range(1,101)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_batch', '{keele_tissue}.{chunk}.txt.gz'), keele_tissue=keele_tissues, chunk=range(1,101)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.{chunk}.txt.gz'), keele_tissue=keele_tissues, chunk=range(1,101)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_transcriptome_variability', '{keele_tissue}.{chunk}.txt.gz'), keele_tissue=keele_tissues, chunk=range(1,101)),
        # Do eqtl as in pape
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_batch', '{keele_tissue}.txt.gz'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_batch', '{keele_tissue}.top_pvalues_per_gene.txt'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_no_covariate', '{keele_tissue}.txt.gz'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_no_covariate', '{keele_tissue}.top_pvalues_per_gene.txt'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability', '{keele_tissue}.txt.gz'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability', '{keele_tissue}.top_pvalues_per_gene.txt'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.txt.gz'), keele_tissue=keele_tissues),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.top_pvalues_per_gene.txt'), keele_tissue=keele_tissues)
    
    
rule get_transcriptome_diversity:
    input: 
        lambda wc: sample_table.loc[wc.dataset, 'path']
    output: 
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', '{dataset}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
        
rule get_cors_with_PCs:
    input: 
        lambda wc: sample_table.loc[wc.dataset, 'path']
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity', '{dataset}', 'cors.txt')
    shell:
        'cd R/; Rscript get_expression_explained_transcriptome_diversity.R {input} {output}'
        
rule expression_association_transcriptome_diversity:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', '{dataset}.txt'),
        lambda wc: sample_table.loc[wc.dataset, 'path']
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations', '{dataset}.txt')
    shell:
        '''
        cd R; Rscript expression_vs_transcriptome_diversity_general.R {input} {output}
        '''
        
#---------------------------------#
# Differential gene expression analysis

rule dge_satantopoulou:
    input:
        os.path.join(config['projectDir'], 'expression_matrices', 'sarantopoulou_{platform}.txt'),
        os.path.join(config['projectDir'], 'expression_matrices', 'sarantopoulou_{platform}_cov.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'sarantopoulou_{platform}_tpm.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'sarantopoulou_{platform}_tmm.txt')
    params:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'sarantopoulou_{platform}')
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'sarantopoulou_{platform}', 'edgeR_results_ILB.Rds')
    shell:
        '''
        cd R
        Rscript differential_expression_sarantopoulou.R {input} {params}
        '''
        
rule dge_nagarajan:
    input:
        os.path.join(config['projectDir'], 'expression_matrices', 'nagarajan.txt'),
        os.path.join(config['projectDir'], 'expression_matrices', 'nagarajan_cov.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'nagarajan_tpm.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'nagarajan_tmm.txt')
    params:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan')
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_genotype.Rds'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_fulvestrant_wt.Rds'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_tamoxifen_wt.Rds'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'differential_gene_expression', 'nagarajan', 'edgeR_results_jq1_wt.Rds')
    shell:
        '''
        cd R
        Rscript differential_expression_nagarajan.R {input} {params}
        '''
 
#---------------------------------#
# Analyzing multiple processing pipelines Arora
# 

rule perform_PCA_arora:
    input:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files'])
    params:
        matrix_dir=os.path.join(config['projectDir'], 'expression_matrices/'),
        project='{project}',
        count_type='{count_type}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_arora', '{project}_{count_type}.Rds'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_arora', '{project}_{count_type}.pca.txt')
    shell:
        '''
        cd R; Rscript PCA_arora.R {input} {params} {output} 
        '''
    

#---------------------------------#
# eqtl analysis keele with GTEx guidelines
# 

rule get_transcriptome_diversity_as_covariates:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'keele_tpm_{keele_tissue}.txt'),
        expression_bed=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}.bed.gz')
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt')
    shell:
        'cd R/; Rscript get_transcriptome_diversity_as_covariates_general.R {input} {output}'

rule merge_non_peer_covs_with_transcriptome_diversity:
    input: 
        trans=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt'),
        peer=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}_covariates.txt')
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates_plus_batch', 'keele_{keele_tissue}.txt')
    shell:
        'Rscript R/rbind.R {output} {input}'
        
rule perform_eqtl_with_no_covariates:
    input: 
        expression_bed=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}.bed.gz'),
        vcf=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db', 'all_keele_mm9.vcf.gz')
    params:
        chunk='{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_no_covariate', '{keele_tissue}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 100 --vcf {input.vcf} --bed {input.expression_bed} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_batch_covariates:
    input: 
        expression_bed=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}.bed.gz'),
        vcf=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db', 'all_keele_mm9.vcf.gz'),
        cov=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}_covariates.txt')
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_batch', '{keele_tissue}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 100 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.cov} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_trasncriptome_diversity_covariates:
    input: 
        expression_bed=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}.bed.gz'),
        vcf=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db', 'all_keele_mm9.vcf.gz'),
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt')
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_transcriptome_variability', '{keele_tissue}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 100 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.cov} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_trasncriptome_diversity_and_batch_covariates:
    input: 
        expression_bed=os.path.join(config['projectDir'], 'expression_matrices', 'keele_{keele_tissue}.bed.gz'),
        vcf=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db', 'all_keele_mm9.vcf.gz'),
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates_plus_batch', 'keele_{keele_tissue}.txt')
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 100 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.cov} --permute 1000 10000 --out {output}
        '''
        
#---------------------------------#
# eqtl analysis keele with as done in their paper
# 

rule perform_eqtl_as_paper_batch_covariates:
    input: 
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt'),
        expression_original=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        gene_info=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        genome_cache_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1/'),
        var_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db/')
    params:
        covariate_mode='batch'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_batch', '{keele_tissue}.txt.gz'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_batch', '{keele_tissue}.top_pvalues_per_gene.txt')
    shell:
        '''
        cd R
        Rscript eqtl_keele_method.R {input} {params} {output}
        '''
        
rule perform_eqtl_as_paper_trasncriptome_diversity_batch_covariates:
    input: 
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt'),
        expression_original=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        gene_info=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        genome_cache_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1/'),
        var_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db/')
    params:
        covariate_mode='transcriptome_diversity+batch'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.txt.gz'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability_and_batch', '{keele_tissue}.top_pvalues_per_gene.txt')
    shell:
        '''
        cd R
        Rscript eqtl_keele_method.R {input} {params} {output}
        '''
        
rule perform_eqtl_as_paper_trasncriptome_diversity_covariates:
    input: 
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt'),
        expression_original=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        gene_info=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        genome_cache_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1/'),
        var_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db/')
    params:
        covariate_mode='transcriptome_diversity'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability', '{keele_tissue}.txt.gz'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_covariate_transcriptome_variability', '{keele_tissue}.top_pvalues_per_gene.txt')
    shell:
        '''
        cd R
        Rscript eqtl_keele_method.R {input} {params} {output}
        '''
rule perform_eqtl_as_paper_no_covariates:
    input: 
        cov=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity_as_covariates', 'keele_{keele_tissue}.txt'),
        expression_original=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        gene_info=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        genome_cache_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1/'),
        var_dir=os.path.join(config['projectDir'], 'expression_datasets', 'keele', 'data_all', 'isvdb_var_db/')
    params:
        covariate_mode='none'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_no_covariate', '{keele_tissue}.txt.gz'),
        os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'keele_paper_eqtls_no_covariate', '{keele_tissue}.top_pvalues_per_gene.txt')
    shell:
        '''
        cd R
        Rscript eqtl_keele_method.R {input} {params} {output}
        '''
