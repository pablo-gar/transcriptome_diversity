# To make expression tables first run:
# cd R; Rscript create_sample_table_DNA_samples.R

import pandas as pd
import os

configfile: '../config.json'
localrules: all, get_chromosome_arms, merge_count_reads_chromosome_arms


#---------------------------------#
# Load tissue sample table
sample_table = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['simple_caller_sample_table_expression']))
col_names = list(sample_table.columns)
col_names[:5] = ['sample_id', 'group', 'path', 'gtex_id_long', 'gtex_id']
sample_table.columns = col_names
sample_table = sample_table.set_index(['sample_id', 'group'])

all_tissues = set(sample_table.index.get_level_values('group'))

# Remove tissues not considered by gtex
all_tissues.remove('Bladder')
all_tissues.remove('Fallopian_Tube')
all_tissues.remove('Cervix_Ectocervix')
all_tissues.remove('Cervix_Endocervix')
#all_tissues=['Whole_Blood']
#---------------------------------#
# Load expression matrix table for fastqtl
expression_table_fastqtl = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_fastqtl']))
expression_table_fastqtl = expression_table_fastqtl.set_index('group')

#---------------------------------#
# Pipeline
# 

rule all:
    input: 
        # Get transcriptome diversity
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt'), group=list(all_tissues)),
        #os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_all_tissues.txt'),
        ## Get expression associations
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'significant.txt.gz'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'nonSignificant.txt.gz'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'GO_results.txt'), group=list(all_tissues)),
        ## Get expression associations no covariates
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'significant.txt.gz'), group=list(all_tissues)),
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'nonSignificant.txt.gz'), group=list(all_tissues)),
        expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'GO_results.txt'), group=list(all_tissues))
        ## Get correlations: transcriptome diversity vs peers
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'significant_cors.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'stats.txt'), group=list(all_tissues)),
        ## Get correlations: transcriptome diversity vs peers
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'significant_cors.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'stats.txt'), group=list(all_tissues)),
        ## Get correlations: transcriptome diversity vs PCs TPM expression
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'PC_matrix.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'PC_stats.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'cors.txt'), group=list(all_tissues)),
        ## Get correlations: transcriptome diversity vs PCs TMM expression
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'PC_matrix.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'PC_stats.txt'), group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'cors.txt'), group=list(all_tissues)),
        # Get variance explained by PEER
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'variance_explained_by_PEER', '{group}.txt'), group=list(all_tissues)),
        #-----------------------------------------------
        # eQTL analysis
        # Dowload gtex data
        #os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrixV8']),
        #os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['sample_annotationV8']),
        ## Get gtex expression matrices by tissue
        #expand(os.path.join(config['projectDir'], 'expression_matrices', 'gtex_{group}.txt'), group=list(all_tissues)),
        # Get transcriptome_diversity as covariates for fastqtl
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_all', '{group}.txt'),  group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates', '{group}.txt'),  group=list(all_tissues)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates_plus_non_peer', '{group}.txt'),  group=list(all_tissues)),
        # Get eQTLs with transcriptome diversity as covariate
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2))
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,201)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,201))
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_PEER_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_nonPEER_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_full_covariate', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate_plus_non_peer', '{group}.{chunk}.txt.gz'), group=list(all_tissues), chunk=range(1,2001))
    
    
rule get_transcriptome_diversity:
    input: 
        exp_mat_file =  os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrixV8']),
    params:
        # Gets sample ids 
        samples = lambda wc: ",".join(list(sample_table.xs(wc.group, level='group')['gtex_id_long'])),
        do_splicing = 'FALSE'
    output: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity.R {input} {params} {output}'
        
rule merge_transcriptome_diversity:
    input: expand(os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt'), group=list(all_tissues))
    output: os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_all_tissues.txt')
    shell: 'cd R; Rscript join_tables.R {output} {input}'
    

rule expression_association_transcriptome_diversity:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        pvalue='0.001',
        do_exon_junctions='FALSE',
        do_exon_junctions_leaf_cutter='FALSE',
        do_covariates='TRUE',
        tissue='{group}'
    output: 
        significant = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'significant.txt.gz'),
        non_significant = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'nonSignificant.txt.gz'),
        go = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations', '{group}', 'GO_results.txt')
    shell: 
        'cd R; Rscript expression_vs_transcriptome_diversity.R {params} {input} {output}'
        
rule expression_association_transcriptome_diversity_no_covariates:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        pvalue='0.001',
        do_exon_junctions='FALSE',
        do_exon_junctions_leaf_cutter='FALSE',
        do_covariates='FALSE',
        tissue='{group}'
    output: 
        significant = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'significant.txt.gz'),
        non_significant = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'nonSignificant.txt.gz'),
        go = os.path.join(config['projectDir'], 'transcriptome_diversity', 'expression_associations_no_covariates', '{group}', 'GO_results.txt')
    shell: 
        'cd R; Rscript expression_vs_transcriptome_diversity.R {params} {input} {output}'
    
rule get_cors_with_PEER:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        tissue='{group}',
        pvalue='0.05',
        gtex_version='v8'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'significant_cors.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'peer_associated_with_transcriptome_diversity', '{group}', 'stats.txt')
    shell:
        'cd R/; Rscript PEER_vs_transcriptome_diversity.R {input} {params} {output}'
        
rule get_cors_with_PCs:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        gtex_version='v8'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity', '{group}', 'cors.txt')
    shell:
        'cd R/; Rscript expression_explained_by_transcriptome_diversity.R {input} {params} {output}'
        
rule get_cors_with_PCs_tmm:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        tissue='{group}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'PCA_var_explained_by_transcriptome_diversity_tmm', '{group}', 'cors.txt')
    shell:
        'cd R/; Rscript get_expression_normalized_explained_transcriptome_diversity.R {input} {params} {output}'

rule get_var_explained_transcriptome_diversity_by_PEER:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity', '{group}.txt')
    params:
        tissue='{group}',
        gtex_version='v8'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'variance_explained_by_PEER', '{group}.txt'),
    shell:
        'cd R/; Rscript PEER_regression_on_transcriptome_diversity.R {input} {params} {output}'
        
#---------------------------------#
# eqtl analysis
# 

rule get_files_from_gtex:
    output:
        exp_mat_file =  os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrixV8']),
        sample_annotation_file = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['sample_annotationV8'])
    shell:
        '''
        wget -O {output.exp_mat_file}.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
        wget -O {output.sample_annotation_file} https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
        
        gunzip {output.exp_mat_file}.gz
        '''

rule get_expression_matrix_by_tissue:
    input:
        exp_mat_file =  os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrixV8']),
        sample_annotation_file = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['sample_annotationV8'])
    params:
        outdir = os.path.join(config['projectDir'], 'expression_matrices')
    output:
        expand(os.path.join(config['projectDir'], 'expression_matrices', 'gtex_{group}.txt'), group=list(all_tissues))
    shell:
        'cd R/; Rscript get_expression_matrices_by_tissue_gtex.R {input} {params}'
        
rule get_transcriptome_diversity_general:
    input: 
        os.path.join(config['projectDir'], 'expression_matrices', 'gtex_{group}.txt')
    output: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_all', '{group}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
 
 
rule get_transcriptome_diversity_as_covariates:
    input: 
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_all', '{group}.txt'),
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path']
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates', '{group}.txt')
    shell:
        'cd R/; Rscript get_transcriptome_diversity_as_covariates.R {input} {output}'

rule merge_non_peer_covs_with_transcriptome_diversity:
    input: 
        trans_diversity = os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates', '{group}.txt'),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates_plus_non_peer', '{group}.txt')
    shell:
        'Rscript R/rbind.R {output} {input}'
        
rule perform_eqtl_with_covariates:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        covariate_file = os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates', '{group}.txt'),
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        #../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --region chr11:5225000-5239000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_covariates_plus:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        covariate_file = os.path.join(config['projectDir'], 'transcriptome_diversity', 'transcriptome_diversity_as_covariates_plus_non_peer', '{group}.txt'),
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_covariate_plus_non_peer', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        #../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --region chr11:5225000-5239000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_no_covariates:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --permute 1000 10000 --out {output}
        #../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --region chr11:5225000-5239000 --vcf {input.vcf} --bed {input.expression_bed} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_full_covariates:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_full_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer} --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_PEER_covariates:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_PEER_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_peer" ]
        then 
            grep -e ID -e InferredCov {input.peer} > {input.peer}.only_peer
        fi
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer}.only_peer --permute 1000 10000 --out {output}
        '''
        
rule perform_eqtl_with_nonPEER_covariates:
    input: 
        expression_bed = lambda wc: expression_table_fastqtl.loc[wc.group, 'path'],
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], 'transcriptome_diversity', 'eqtls_nonPEER_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_non_peer" ]
        then 
            grep -v InferredCov {input.peer} > {input.peer}.only_non_peer
        fi
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer}.only_non_peer --permute 1000 10000 --out {output}
        '''
