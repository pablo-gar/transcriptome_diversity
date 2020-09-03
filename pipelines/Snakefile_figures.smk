# To make expression tables first run:
# cd R; Rscript create_sample_table_DNA_samples.R

import pandas as pd
import os

configfile: '../config.json'
localrules: all

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
        # Get lin data
        'figures/data/lin_tpm.txt',
        'figures/data/lin_tpm.txt',
        'figures/data/lin_covariates.txt',
        'figures/data/lin_tpm_transcriptome_diversity.txt',
        'figures/data/lin_tmm_transcriptome_diversity.txt',
        'figures/data/lin_expression_associations_tpm.txt',
        'figures/data/lin_expression_associations_tmm.txt',
        'figures/data/lin_tpm_pca/cors.txt',
        'figures/data/lin_tpm_pca/PC_matrix.txt',
        'figures/data/lin_tpm_pca/PC_stats.txt',
        'figures/data/lin_tmm_pca/cors.txt',
        'figures/data/lin_tmm_pca/PC_matrix.txt',
        'figures/data/lin_tmm_pca/PC_stats.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechPico_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_Illumina_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechV4_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechPico_tpm_transcriptome_diversity.txt',
        'figures/data/sarantopoulou/sarantopoulou_Illumina_tpm_transcriptome_diversity.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechV4_tpm_transcriptome_diversity.txt'
        
rule get_lin_data:
    input:
        tpm_exp=sample_table.loc['lin_tpm', 'path'],
        tmm_exp=sample_table.loc['lin_tmm', 'path'],
        tpm_trans_diversity=expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'transcriptome_diversity', 'lin_tpm.txt')),
        covariates=os.path.join(config['projectDir'], 'expression_matrices', 'lin_covariates.txt')
    output:
        tpm_exp='figures/data/lin_tpm.txt',
        tmm_exp='figures/data/lin_tmm.txt',
        tpm_trans_diversity='figures/data/lin_tpm_transcriptome_diversity.txt',
        covariates='figures/data/lin_covariates.txt'
    shell:
        '''
        cp {input.tpm_exp} {output.tpm_exp}
        cp {input.tmm_exp} {output.tmm_exp}
        cp {input.tpm_trans_diversity} {output.tpm_trans_diversity}
        cp {input.covariates} {output.covariates}
        '''
        
rule get_lin_expression_associations:
    input:
        tpm_associations=expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations_from_raw', 'expression_associations_tpm', 'lin.txt')),
        tmm_associations=expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations_from_raw', 'expression_associations_tmm', 'lin.txt'))
    output:
        tpm_associations='figures/data/lin_expression_associations_tpm.txt',
        tmm_associations='figures/data/lin_expression_associations_tmm.txt'
    shell:
        '''
        cp {input.tpm_associations} {output.tpm_associations}
        cp {input.tmm_associations} {output.tmm_associations}
        '''
        
rule get_lin_PCA_results:
    input:
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity_from_raw', '{exp_type}', 'lin', 'PC_matrix.txt'), exp_type=['tpm', 'tmm']),
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity_from_raw', '{exp_type}', 'lin', 'PC_stats.txt'), exp_type=['tpm', 'tmm']),
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity_from_raw', '{exp_type}', 'lin', 'cors.txt'), exp_type=['tpm', 'tmm'])
    params:
        tpm_in=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity_from_raw', 'tpm', 'lin'),
        tmm_in=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'PCA_var_explained_by_transcriptome_diversity_from_raw', 'tmm', 'lin'),
        tpm_out='figures/data/lin_tpm_pca/',
        tmm_out='figures/data/lin_tmm_pca/'
    output:
        'figures/data/lin_tpm_pca/cors.txt',
        'figures/data/lin_tpm_pca/PC_matrix.txt',
        'figures/data/lin_tpm_pca/PC_stats.txt',
        'figures/data/lin_tmm_pca/cors.txt',
        'figures/data/lin_tmm_pca/PC_matrix.txt',
        'figures/data/lin_tmm_pca/PC_stats.txt'
    shell:
        '''
        cp {params.tpm_in}/* {params.tpm_out}
        cp {params.tmm_in}/* {params.tmm_out}
        '''
        
rule get_sarantopoulou_data:
    input:
        pico_cov=os.path.join(config['projectDir'], 'expression_matrices', 'sarantopoulou_pico_cov.txt'),
        truseq_cov=os.path.join(config['projectDir'], 'expression_matrices', 'sarantopoulou_truseq_cov.txt'),
        v4_cov=os.path.join(config['projectDir'], 'expression_matrices', 'sarantopoulou_v4_cov.txt'),
        pico_trans=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations_from_raw', 'transcriptome_diversity', 'sarantopoulou_pico.txt'),
        truseq_trans=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations_from_raw', 'transcriptome_diversity', 'sarantopoulou_truseq.txt'),
        v4_trans=os.path.join(config['projectDir'], 'transcriptome_diversity_other_datasets', 'expression_associations_from_raw', 'transcriptome_diversity', 'sarantopoulou_v4.txt')
    output:
        pico_cov='figures/data/sarantopoulou/sarantopoulou_ClontechPico_covariates.txt',
        truseq_cov='figures/data/sarantopoulou/sarantopoulou_Illumina_covariates.txt',
        v4_cov='figures/data/sarantopoulou/sarantopoulou_ClontechV4_covariates.txt',
        pico_trans='figures/data/sarantopoulou/sarantopoulou_ClontechPico_tpm_transcriptome_diversity.txt',
        truseq_trans='figures/data/sarantopoulou/sarantopoulou_Illumina_tpm_transcriptome_diversity.txt',
        v4_trans='figures/data/sarantopoulou/sarantopoulou_ClontechV4_tpm_transcriptome_diversity.txt'
    shell:
        '''
        cp {input.pico_cov} {output.pico_cov}
        cp {input.truseq_cov} {output.truseq_cov}
        cp {input.v4_cov} {output.v4_cov}
        cp {input.pico_trans} {output.pico_trans}
        cp {input.truseq_trans} {output.truseq_trans}
        cp {input.v4_trans} {output.v4_trans}
        '''
        
        
        
        
