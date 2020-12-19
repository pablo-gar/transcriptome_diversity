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
sample_table_original = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files']))
sample_table_original = sample_table_original.set_index('id')

sample_table_raw = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_for_expression_associations']))
sample_table_raw = sample_table_raw.set_index('id')


#---------------------------------#
# Pipeline
# 

rule all:
    input: 
        # Figure 2 - simple example
        'figures/figures/figure_2_example/association_example_tpm.pdf',
        'figures/figures/figure_2_example/association_example_tmm.pdf',
        #
        'figures/figures/figure_2_example/all_expression_associations_lin_N.pdf',
        'figures/figures/figure_2_example/all_expression_associations_lin_directionality_tpm.pdf',
        'figures/figures/figure_2_example/all_expression_associations_lin_directionality_tmm.pdf',
        #
        'figures/figures/figure_2_example/lin_pcs_PC1_PC2_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_PC1_PC2_tmm.pdf',
        'figures/figures/figure_2_example/lin_pcs_var_explained_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_var_explained_tmm.pdf',
        'figures/figures/figure_2_example/lin_pcs_cor_trans_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_cor_trans_tmm.pdf',
        #
        'figures/figures/figure_2_example/lin_peer_cors.pdf',
        # Figure 3 - associations with metadata/artifacts
        #'figures/figures/figure_3_features_associated/lin_replicates_strains.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_replicates.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_fly_number.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_sex.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_read_depth.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_model_pvalues.pdf',
        #'figures/figures/figure_3_features_associated/lin_replicates_final_combined_features.pdf',
        ##
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_individual.pdf',
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_untreated.pdf',
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_treat_ILB.pdf',
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_both_treat_untreat.pdf',
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_read_depth.pdf',
        #'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_read_depth_rankit.pdf',
        ## Figure 4 - big consortia results (GTEx, TCGA)
        #'figures/figures/figure_4_big_consortia/gtex_PC_cors_tpm.pdf',
        #'figures/figures/figure_4_big_consortia/gtex_PC_cors_tmm.pdf',
        ##
        #'figures/figures/figure_4_big_consortia/gtex_expression_associations_tmm.pdf',
        #'figures/figures/figure_4_big_consortia/gtex_expression_associations_tpm.pdf',
        #expand(os.path.join('figures/data/keele', 'eqtl_result', '{covariate_mode}', '{keele_tissue}.txt.gz'), covariate_mode=['transcriptome_diversity+batch', 'batch'], keele_tissue=keele_tissues),
        ## Figure 4.1 - GTEx peer
        'figures/figures/figure_4.1_GTEx_PEER/gtex_PEER_cors.pdf',
        
        
rule get_all_data:
    input:
        'figures/data/lin/lin_tpm.txt',
        #'figures/data/lin/lin_tpm.txt',
        #'figures/data/lin/lin_covariates.txt',
        #'figures/data/lin/lin_tpm_transcriptome_diversity.txt',
        #'figures/data/lin/lin_expression_associations_tpm.txt',
        #'figures/data/lin/lin_expression_associations_tmm.txt',
        #'figures/data/lin/lin_tpm_pca/cors.txt',
        #'figures/data/lin/lin_tpm_pca/PC_matrix.txt',
        #'figures/data/lin/lin_tpm_pca/PC_stats.txt',
        #'figures/data/lin/lin_tmm_pca/cors.txt',
        #'figures/data/lin/lin_tmm_pca/PC_matrix.txt',
        #'figures/data/lin/lin_tmm_pca/PC_stats.txt',
        #'figures/data/lin/lin_peer_cors.txt',
        #'figures/data/sarantopoulou/sarantopoulou_ClontechPico_covariates.txt',
        #'figures/data/sarantopoulou/sarantopoulou_Illumina_covariates.txt',
        #'figures/data/sarantopoulou/sarantopoulou_ClontechV4_covariates.txt',
        #'figures/data/sarantopoulou/sarantopoulou_ClontechPico_tpm_transcriptome_diversity.txt',
        #'figures/data/sarantopoulou/sarantopoulou_Illumina_tpm_transcriptome_diversity.txt',
        #'figures/data/sarantopoulou/sarantopoulou_ClontechV4_tpm_transcriptome_diversity.txt',
        #'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tpm/gtex_Whole_Blood/cors.txt',
        #'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tmm/gtex_Whole_Blood/cors.txt',
        #'figures/data/gtex/expression_associations/tpm/gtex_Whole_Blood.txt',
        #'figures/data/gtex/expression_associations/tmm/gtex_Whole_Blood.txt',
        #'figures/data/gtex/eqtl/eqtl_results.txt',
        
#---------------------------------------------------------------
# Make figures


# Fig 2 - Simple example
 
rule make_lin_example:
    input:
        'figures/data/lin/lin_tpm.txt',
        'figures/data/lin/lin_tpm.txt',
        'figures/data/lin/lin_tpm_transcriptome_diversity.txt'
    output:
        'figures/figures/figure_2_example/association_example_tpm.pdf',
        'figures/figures/figure_2_example/association_example_tmm.pdf'
    shell:
        'cd figures/scripts/; Rscript 1_plot_example_lin.R'
        
rule make_lin_associations_with_expression_simple:
    input:
        'figures/data/lin/lin_expression_associations_tpm.txt',
        'figures/data/lin/lin_expression_associations_tmm.txt'
    output:
        'figures/figures/figure_2_example/all_expression_associations_lin_N.pdf',
        'figures/figures/figure_2_example/all_expression_associations_lin_directionality_tpm.pdf',
        'figures/figures/figure_2_example/all_expression_associations_lin_directionality_tmm.pdf'
    shell:
        'cd figures/scripts/; Rscript 1_plot_expression_associations_lin_general_numbers.R'
        
rule make_lin_PCA:
    input:
        'figures/data/lin/lin_covariates.txt',
        'figures/data/lin/lin_tpm_pca/cors.txt',
        'figures/data/lin/lin_tpm_pca/PC_matrix.txt',
        'figures/data/lin/lin_tpm_pca/PC_stats.txt',
        'figures/data/lin/lin_tmm_pca/cors.txt',
        'figures/data/lin/lin_tmm_pca/PC_matrix.txt',
        'figures/data/lin/lin_tmm_pca/PC_stats.txt'
    output:
        'figures/figures/figure_2_example/lin_pcs_PC1_PC2_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_PC1_PC2_tmm.pdf',
        'figures/figures/figure_2_example/lin_pcs_var_explained_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_var_explained_tmm.pdf',
        'figures/figures/figure_2_example/lin_pcs_cor_trans_tpm.pdf',
        'figures/figures/figure_2_example/lin_pcs_cor_trans_tmm.pdf'
    shell:
        'cd figures/scripts/; Rscript 1_plot_pca_lin.R'
        
rule make_lin_peer_cors:
    input:
        'figures/data/lin/lin_peer_cors.txt'
    output:
        'figures/figures/figure_2_example/lin_peer_cors.pdf'
    shell:
        'cd figures/scripts/; Rscript 1_lin_PEER.R'
        
        
        
# Fig 3 - technical and biological artifacts associated with transcriptome diversity
        
rule make_lin_artifacts:
    input:
        'figures/data/lin/lin_covariates.txt',
        'figures/data/lin/lin_tpm_transcriptome_diversity.txt'
    output:
        'figures/figures/figure_3_features_associated/lin_replicates_strains.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_replicates.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_fly_number.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_sex.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_read_depth.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_model_pvalues.pdf',
        'figures/figures/figure_3_features_associated/lin_replicates_final_combined_features.pdf'
    shell:
        'cd figures/scripts/; Rscript 2_lin_artifacts.R'
        
rule make_sarantopoulou_platform_differences:
    input:
        'figures/data/sarantopoulou/sarantopoulou_ClontechPico_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_Illumina_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechV4_covariates.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechPico_tpm_transcriptome_diversity.txt',
        'figures/data/sarantopoulou/sarantopoulou_Illumina_tpm_transcriptome_diversity.txt',
        'figures/data/sarantopoulou/sarantopoulou_ClontechV4_tpm_transcriptome_diversity.txt'
    output:
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_individual.pdf',
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_untreated.pdf',
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_treat_ILB.pdf',
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_Platform_both_treat_untreat.pdf',
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_read_depth.pdf',
        'figures/figures/figure_3_features_associated/2_sarantopoulou_platforms_read_depth_rankit.pdf'
    shell:
        'cd figures/scripts/; Rscript 2_sarantopoulou_plaforms_info.R'
        
        
# Fig 4 - big consortia figures (GTEx, TCGA)
 
# Correlations of PCs with transcriptome diversity
rule make_gtex_PC_cors_transcriptome_diversity:
    input:
        'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tpm/gtex_Whole_Blood/cors.txt',
        'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tmm/gtex_Whole_Blood/cors.txt'
    output:
        'figures/figures/figure_4_big_consortia/gtex_PC_cors_tpm.pdf',
        'figures/figures/figure_4_big_consortia/gtex_PC_cors_tmm.pdf'
    shell:
        'cd figures/scripts/; Rscript 3_gtex_PCs_cors_trans.R'
        
        
# Number of genes whose expression is associated with transcriptome diversity
rule make_expression_associations_gtex:
    input:
        'figures/data/gtex/expression_associations/tpm/gtex_Whole_Blood.txt',
        'figures/data/gtex/expression_associations/tmm/gtex_Whole_Blood.txt'
    output:
        'figures/figures/figure_4_big_consortia/gtex_expression_associations_tmm.pdf',
        'figures/figures/figure_4_big_consortia/gtex_expression_associations_tpm.pdf'
    shell:
        'cd figures/scripts/; Rscript 3_gtex_expression_associations.R'
        
# Fig 4.1 - GTEx PEER
 
# Correlations of PCs with transcriptome diversity
rule make_GTEx_peer:
    input:
        'figures/data/gtex/PEER/Whole_Blood/PEER_transcriptome_diversity.txt'
    output:
        'figures/figures/figure_4.1_GTEx_PEER/gtex_PEER_cors.pdf',
    shell:
        'cd figures/scripts/; Rscript 4_gtex_PEER_cors_trans_diversity.R'
        
        
        
#---------------------------------------------------------------
# Transfer data
        
rule get_lin_data:
    input:
        #tpm_exp=os.path.join(config['projectDir'], sample_table.loc['lin_tpm', 'path']),
        #tmm_exp=os.path.join(config['projectDir'], sample_table.loc['lin_tmm', 'path']),
        tpm_exp=os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tpm', 'lin.txt'),
        tmm_exp=os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tmm', 'lin.txt'),
        tpm_trans_diversity=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'lin.txt'),
        peer_cors=os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', 'lin', 'cors.txt'),
        covariates=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_covariates.txt')
    output:
        tpm_exp='figures/data/lin/lin_tpm.txt',
        tmm_exp='figures/data/lin/lin_tmm.txt',
        tpm_trans_diversity='figures/data/lin/lin_tpm_transcriptome_diversity.txt',
        peer_cors='figures/data/lin/lin_peer_cors.txt',
        covariates='figures/data/lin/lin_covariates.txt'
    shell:
        '''
        cp {input.tpm_exp} {output.tpm_exp}
        cp {input.tmm_exp} {output.tmm_exp}
        cp {input.tpm_trans_diversity} {output.tpm_trans_diversity}
        cp {input.covariates} {output.covariates}
        cp {input.peer_cors} {output.peer_cors}
        '''
        
rule get_lin_expression_associations:
    input:
        tpm_associations=os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', 'tpm', 'lin.txt'),
        tmm_associations=os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', 'tmm', 'lin.txt')
    output:
        tpm_associations='figures/data/lin/lin_expression_associations_tpm.txt',
        tmm_associations='figures/data/lin/lin_expression_associations_tmm.txt'
    shell:
        '''
        cp {input.tpm_associations} {output.tpm_associations}
        cp {input.tmm_associations} {output.tmm_associations}
        '''
        
rule get_lin_PCA_results:
    input:
        expand(os.path.join(config['projectDir'], config['expression_associations_dir'], 'PCA_var_explained_by_transcriptome_diversity', '{exp_type}', 'lin', 'PC_matrix.txt'), exp_type=['tpm', 'tmm']),
        expand(os.path.join(config['projectDir'], config['expression_associations_dir'], 'PCA_var_explained_by_transcriptome_diversity', '{exp_type}', 'lin', 'PC_stats.txt'), exp_type=['tpm', 'tmm']),
        expand(os.path.join(config['projectDir'], config['expression_associations_dir'], 'PCA_var_explained_by_transcriptome_diversity', '{exp_type}', 'lin', 'cors.txt'), exp_type=['tpm', 'tmm'])
    params:
        tpm_in=os.path.join(config['projectDir'], config['expression_associations_dir'], 'PCA_var_explained_by_transcriptome_diversity', 'tpm', 'lin'),
        tmm_in=os.path.join(config['projectDir'], config['expression_associations_dir'], 'PCA_var_explained_by_transcriptome_diversity', 'tmm', 'lin'),
        tpm_out='figures/data/lin/lin_tpm_pca/',
        tmm_out='figures/data/lin/lin_tmm_pca/'
    output:
        'figures/data/lin/lin_tpm_pca/cors.txt',
        'figures/data/lin/lin_tpm_pca/PC_matrix.txt',
        'figures/data/lin/lin_tpm_pca/PC_stats.txt',
        'figures/data/lin/lin_tmm_pca/cors.txt',
        'figures/data/lin/lin_tmm_pca/PC_matrix.txt',
        'figures/data/lin/lin_tmm_pca/PC_stats.txt'
    shell:
        '''
        cp {params.tpm_in}/* {params.tpm_out}
        cp {params.tmm_in}/* {params.tmm_out}
        '''
        
rule get_sarantopoulou_data:
    input:
        pico_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'],  'sarantopoulou_pico_cov.txt'),
        truseq_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'],  'sarantopoulou_truseq_cov.txt'),
        v4_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'],  'sarantopoulou_v4_cov.txt'),
        pico_trans=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'sarantopoulou_pico.txt'),
        truseq_trans=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'sarantopoulou_truseq.txt'),
        v4_trans=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'sarantopoulou_v4.txt')
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

rule get_gtex_PC_data:
    # Only specify whole blood, but it will get everything
    input:
        os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw','tpm', 'gtex_Whole_Blood', 'cors.txt'),
        os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw','tmm', 'gtex_Whole_Blood', 'cors.txt')
    params:
        in_folder=os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw'),
        out='figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/'
    output:
        'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tpm/gtex_Whole_Blood/cors.txt',
        'figures/data/gtex/PCA_var_explained_by_transcriptome_diversity/tmm/gtex_Whole_Blood/cors.txt'
    shell:
        '''
        cp -r {params.in_folder}/tpm/gtex_* {params.out}/tpm/
        cp -r {params.in_folder}/tmm/gtex_* {params.out}/tmm/
        '''
        
rule get_gtex_expression_associations:
    input:
        tpm_associations=os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', 'tpm', 'gtex_Whole_Blood.txt'),
        tmm_associations=os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', 'tmm', 'gtex_Whole_Blood.txt')
    output:
        tpm_associations='figures/data/gtex/expression_associations/tpm/gtex_Whole_Blood.txt',
        tmm_associations='figures/data/gtex/expression_associations/tmm/gtex_Whole_Blood.txt'
    shell:
        '''
        cp $(dirname {input.tpm_associations})/gtex_* $(dirname {output.tpm_associations})/ 
        cp $(dirname {input.tmm_associations})/gtex_* $(dirname {output.tmm_associations})/ 
        '''
        
rule get_gtex_eqtl_results:
    input:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'all_results', 'gtex_eqtl_results.txt')
    output:
        'figures/data/gtex/eqtl/eqtl_results.txt'
    shell:
        '''
        cp {input} {output}
        '''
    
rule get_gtex_PEER_cors:
    input:
        cors=os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', 'Whole_Blood', 'cors.txt'),
        stats=os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', 'Whole_Blood', 'stats.txt'),
        peer_trans=os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', 'Whole_Blood', 'PEER_transcriptome_diversity.txt')
    output:
        cors='figures/data/gtex/PEER/Whole_Blood/cors.txt',
        stats='figures/data/gtex/PEER/Whole_Blood/stats.txt',
        peer_trans='figures/data/gtex/PEER/Whole_Blood/PEER_transcriptome_diversity.txt'
    shell:
        '''
        cp {input.cors} {output.cors}
        cp {input.stats} {output.stats}
        cp {input.peer_trans} {output.peer_trans}
        '''
        
rule get_keele_eqtl_results:
    input:
        expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'keele', '{covariate_mode}', '{keele_tissue}.txt.gz'), covariate_mode=['transcriptome_diversity+batch', 'batch'], keele_tissue=keele_tissues)
    params:
        in_folder=os.path.join(config['projectDir'], config['eqtl_dir'], 'keele'),
        out_folder=os.path.join('figures/data/keele', 'eqtl_result')
    output:
        expand(os.path.join('figures/data/keele', 'eqtl_result', '{covariate_mode}', '{keele_tissue}.txt.gz'), covariate_mode=['transcriptome_diversity+batch', 'batch'], keele_tissue=keele_tissues)
    shell:
        '''
        cp -r {params.in_folder}/* {params.out_folder}
        '''
