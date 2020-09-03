# To make expression tables first run:
# cd R; Rscript create_sample_table_DNA_samples.R
# 
# ------------------------------------------------------
# Required files (manual downloads)
# 
# Convert this xlsx to tsv (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60314/suppl/GSE60314_GEO_run_summary.xlsx)
# os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin/GSE60314_GEO_run_summary_metadata.txt/')
# 
# Dowload this file from : https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE140872&format=file
# os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'GSE140872_RAW.tar')
# 
# The GTEx v8 metadata (phenotypes) were downloaded from dbagp (they are authorized access only)
# os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'GTEx_Phenotypes_v8/'),

import os

configfile: '../config.json'
localrules: all

keele_tissues=['kidney', 'liver', 'lung']

#---------------------------------#
# Pipeline
# 

rule all:
    input: 
        # Download Arora et al. 2020
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora', 'OriginalTCGAGTExData', 'file_md5sums.txt'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora', 'OriginalTCGAGTExData', 'SE_objects/gtex_mskcc_batch_log2_TPM.RData'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_arora.txt'),
        #Download keele 
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1', 'genotype', 'chr1', 'markers.RData' ),
        expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.rds'), chr=[str(i) for i in range(1,20)] + ['X']),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'GEO_matrices','GSM4189358_CC068_lung_gene_quant.txt.gz'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_paper', 'pgen.1008537.s045.R'),
        #Create keele vcf files
        #expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.rds'), chr=[str(i) for i in range(1,20)] + ['X']),
        #os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', 'all_keele.#vcf.gz'),
        #Create keele standard expression and covariate files
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.txt'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.bed.gz'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}_covariates.txt'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_tpm_{keele_tissue}.txt'), keele_tissue=keele_tissues),
        #Create keele expression look-up table
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_keele.txt'),
        #
        # Download lin
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_5_57_HTSeq_raw_read_counts.txt.gz'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_6_01_HTSeq_raw_read_counts.txt.gz'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_GEO_run_summary.xlsx'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_dme5_57_ERCC.gtf.gz'),
        # Create lin matrices
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_counts.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_tmm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_covariates.txt'),
        # Download sarantopoulou
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_sarantopoulou.txt'),
        # Download Nagarajan data
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_nagarajan.txt'),
        # Dowload gtex
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'gencode.v26.GRCh38.genes.gtf'),
        os.path.join(config['projectDir'], config['expression_datasets_gtex_bed_dir'], 'Whole_Blood.v8.normalized_expression.bed.gz'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_gtex_bed']),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'gtex_Whole_Blood_counts.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'gtex_Whole_Blood.txt'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['sample_annotationV8_GTEx']),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'Whole_Blood.v8.covariates.txt'),
        ## Create look up table for expression association (raw counts)
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_raw_counts.txt'),
        ## Get mm9 id-to-name conversion table
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.id_converstion_table.txt')

##################################################
## Getting reading look up table for raw counts

rule create_lookup_table_for_expression:
    input:
        #Create keele standard expression and covariate files
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.txt'), keele_tissue=keele_tissues),
        # Create lin matrices
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_counts.txt'),
        # Download sarantopoulou
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_sarantopoulou.txt'),
        # Download Nagarajan data
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_nagarajan.txt'),
        # Dowload gtex
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'),
    params:
        root_dir=config['projectDir'],
        exp_dir=config['expression_matrices_dir'],
        aux_dir=config['auxiliary_Files']['dir'],
        datasets_dir=config['expression_datasets_dir']
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_raw_counts.txt')
    shell:
        '''
        Rscript R/create_lookup_table_for_raw_counts.R {params} {output}
        '''
                

##################################################
## GTEx 

rule download_gtex_raw_counts: 
    params:
        os.path.join(config['projectDir'], config['expression_matrices_dir'])
    output:
        exp=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'),
        anno=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'),
        gtf=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'gtex', 'gencode.v26.GRCh38.genes.gtf'),
        blood=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'gtex_Whole_Blood_counts.txt')
    shell:
        '''
        wget -O {output.exp} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
        wget -O {output.anno} https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
        wget -O {output.gtf} https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf
        
        cd R
        Rscript process_gtex_data_raw_counts.R {output.exp} {output.anno} {params}
        '''

rule download_gtex_bed:
    output:
        blood=os.path.join(config['projectDir'], config['expression_datasets_gtex_bed_dir'], 'Whole_Blood.v8.normalized_expression.bed.gz'),
        sample_matrix=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_gtex_bed'])
    params:
        relative_path=config['expression_datasets_gtex_bed_dir']
    shell:
        '''
        # Downlad data
        dir=$(dirname {output.blood})
        wget -O $dir/gtex.tar https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar
        tar -C $dir -xvf $dir/gtex.tar
        mv $dir/GTEx_Analysis_v8_eQTL_expression_matrices/* $dir/
        touch $dir/*tbi
        rm -r $dir/GTEx_Analysis_v8_eQTL_expression_matrices/
        rm $dir/gtex.tar
        
        # Make sample table
        cd  $dir
        for i in *gz
        do
            paste <(perl -pe 's/(.+?)\..+/$1/' <(echo $i)) <(echo {params.relative_path}/$i)
        done > {output.sample_matrix}.temp
        cat <(echo -e 'group\\tpath') {output.sample_matrix}.temp > {output.sample_matrix} && rm {output.sample_matrix}.temp
        '''
        
rule dowload_gtex_tpm_and_annotation:
    params:
        temp=os.path.join(config['projectDir'], config['temp']),
        root_dir=config['projectDir'],
        relative_dir=config['expression_matrices_dir']
    output:
        blood=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'gtex_Whole_Blood.txt'),
        sample_annotation_file=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['sample_annotationV8_GTEx']),
        tissue_table=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_gtex_tpm'])
    shell:
        '''
        wget -O {params.temp}/gtex.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
        wget -O {output.sample_annotation_file} https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
        
        cd R/; Rscript R/process_arora_data_gtex_table.txt {params.temp}/gtex.gz {output.sample_annotation_file} {output.tissue_table} {params.root_dir} {params.relative_dir}
        
        rm {params.temp}/gtex.gz
        '''
        
rule dowload_gtex_covariates:
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'Whole_Blood.v8.covariates.txt'),
    params:
        temp=os.path.join(config['projectDir'], config['temp']),
    shell:
        '''
        wget -O {params.temp}/covariates_gtex.tar.gz https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz
        
        tar -C $(dirname {output}) -xvzf {params.temp}/covariates_gtex.tar.gz
        
        mv $(dirname {output})/GTEx_Analysis_v8_eQTL_covariates/* $(dirname {output})
        
        rm {params.temp}/covariates_gtex.tar.gz
        rm -r $(dirname {output})/GTEx_Analysis_v8_eQTL_covariates/
        '''
        
##################################################
## Arora

rule download_arora:
    params:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora')
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora', 'OriginalTCGAGTExData', 'file_md5sums.txt'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora', 'OriginalTCGAGTExData', 'SE_objects/gtex_mskcc_batch_log2_TPM.RData')
    shell:
        '''
        cd {params}
        
        wget --recursive -nH --cut-dirs=1 https://s3-us-west-2.amazonaws.com/fh-pi-holland-e/OriginalTCGAGTExData/index.html
        '''
        
rule create_arora_matrices:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'arora', 'OriginalTCGAGTExData', 'SE_objects/gtex_mskcc_batch_log2_TPM.RData'),
    params:
        out_dir_base=config['projectDir'],
        out_dir_final=config['expression_matrices_dir']
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_arora.txt')
    shell:
        '''
        mkdir -p {params}
        cd R
        Rscript process_arora_data.R $(dirname {input}) {params} {output}
        '''
        
        
        
        
##################################################
## Keele
        
rule download_keele:
    params:
        root=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele'),
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1', 'genotype', 'chr1', 'markers.RData' ),
        expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.rds'), chr=[str(i) for i in range(1,20)] + ['X']),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_paper', 'pgen.1008537.s045.R')
    shell:
        '''
        mkdir -p {params.root}
        cd {params.root}
        
        mkdir -p data_all
        
        wget -O keele.zip https://ndownloader.figshare.com/articles/9985514/versions/3
        wget -O keele_paper.zip https://ndownloader.figstatic.com/articles/11675130/versions/1
        
        unzip keele.zip -d data_all
        
        unzip data_all/cc_genome_cache_full_l2_0.1.zip -d data_all/
        
        mkdir -p data_all/isvdb_var_db
        tar -C data_all/isvdb_var_db -xvzf data_all/isvdb_var_db.tar.gz 
        
        mkdir -p data_paper
        unzip keele_paper.zip -d data_paper/
        
        rm keele.zip
        rm data_all/cc_genome_cache_full_l2_0.1.zip
        rm data_all/isvdb_var_db.tar.gz 
        rm keele_paper.zip
        '''


# This requires to manually download GSE140872_RAW.tar (see at the top for the link)
rule uncompress_keele_geo:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'GSE140872_RAW.tar')
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'GEO_matrices','GSM4189358_CC068_lung_gene_quant.txt.gz'),
    shell:
        '''
        tar -xf {input} -C $(dirname {output})
        '''
        

rule create_keele_vcfs:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.rds')
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.vcf')
    shell:
        '''
        module load bcftools
        module load htslib
        Rscript R/process_keele_isvdb_to_vcf.R {input} {output}.temp
        bcftools sort -o {output} {output}.temp
        tabix {output}
        rm {output}.temp
        '''
        
rule merge_vcfs:
    input:
        expand(os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', '{chr}.vcf'), chr=[str(i) for i in range(1,20)] + ['X']),
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db', 'all_keele.vcf.gz')
    shell:
        '''
        module load bcftools
        bcftools concat -O z -o {output} {input}
        bcftools index {output}
        '''
        
rule create_keele_matrices:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.id_converstion_table.txt')
    params:
        bed=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.bed'),
    output:
        exp_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.txt'),
        cov=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}_covariates.txt'),
        bed=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.bed.gz'),
    shell:
        '''
        module load htslib
        cd R
        Rscript process_keele_create_expression_and_covariates.R {input} {output.exp_mat} {output.cov} {params}
        bgzip {params.bed}
        tabix {output.bed}
        '''

rule create_keele_matrices_tpm:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'GEO_matrices','GSM4189358_CC068_lung_gene_quant.txt.gz'),
    params:
        prefix='keele_temp_',
        out_dir=os.path.join(config['projectDir'], config['expression_matrices_dir'])
    output: 
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_temp_{aaa}.txt'), aaa=keele_tissues)
    shell:
        '''
        cd R
        Rscript process_keele_create_expression_tpm.R $(dirname {input}) {params}
        '''
        
rule convert_keele_names_to_ids:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_temp_{tissue}.txt'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.id_converstion_table.txt')
    output: 
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_tpm_{tissue}.txt')
    shell:
        '''
        cd R
        Rscript process_keele_convert_names_to_ids.R {input} {output}
        '''
        
rule create_keele_matrix_look_up_table:
    input:
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_tpm_{keele_tissue}.txt'), keele_tissue=keele_tissues),
        expand(os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.txt'), keele_tissue=keele_tissues)
    params:
        ' '.join(expand('keele_tpm_{keele_tissue}', keele_tissue=keele_tissues) + expand('keele_{keele_tissue}', keele_tissue=keele_tissues))
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_keele.txt')
    shell:
        '''
        echo "{input}" | perl -pe 's/ /\\n/g' > {output}.1
        echo "{params}" | perl -pe 's/ /\\n/g' > {output}.2
        
        paste {output}.2 {output}.2 {output}.1 > {output}.3
        
        cat <(echo -e 'id\\tname\\tpath') {output}.3 > {output}
        
        rm {output}.1 {output}.2 {output}.3
        '''
        
        
        
        
##################################################
## Lin

rule download_lin:
    output:
        raw_1=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_5_57_HTSeq_raw_read_counts.txt.gz'),
        raw_2=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_6_01_HTSeq_raw_read_counts.txt.gz'),
        summary=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_GEO_run_summary.xlsx'),
        gtf=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_dme5_57_ERCC.gtf.gz')
        
    shell:
        '''
        wget -O {output.raw_1} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60314/suppl/GSE60314_5_57_HTSeq_raw_read_counts.txt.gz
        wget -O {output.raw_2} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60314/suppl/GSE60314_6_01_HTSeq_raw_read_counts.txt.gz
        wget -O {output.summary} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60314/suppl/GSE60314_GEO_run_summary.xlsx
        wget -O {output.gtf} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60314/suppl/GSE60314_dme5_57_ERCC.gtf.gz
        '''
        
rule create_lin_matrices:
    input: 
        raw_1=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_5_57_HTSeq_raw_read_counts.txt.gz'),
        summary=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_GEO_run_summary_metadata.txt')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_counts.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_covariates.txt')
    shell:
        '''
        cd R; Rscript process_lin_create_expression.R {input} {output}
        '''
        
rule create_lin_matrices_tpm:
    input: 
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_counts.txt'),
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'lin', 'GSE60314_dme5_57_ERCC.gtf.gz')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin.txt')
    shell:
        '''
        cd R; Rscript get_TPM.R {input} {output}
        '''
        
rule create_lin_matrices_tmm:
    input: 
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_counts.txt')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'lin_tmm.txt')
    shell:
        '''
        cd R; Rscript get_tmm.R {input} {output}
        '''



##################################################
## sarantopoulou

rule download_sarantopoulou:
    output:
        v4=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.V4.txt.gz'),
        pico=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Pico.R.txt.gz'),
        truseq=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Truseq.txt.gz')
    shell:
        '''
        wget -O {output.v4} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.V4.txt.gz
        wget -O {output.pico} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Pico.R.txt.gz
        wget -O {output.truseq} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Truseq.txt.gz
        '''
        
rule create_sarantopoulou_matrices_counts:
    input:
        v4=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.V4.txt.gz'),
        pico=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Pico.R.txt.gz'),
        truseq=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'sarantopoulou', 'GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Truseq.txt.gz')
    output:
        v4_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4.txt'),
        v4_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4_cov.txt'),
        pico_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico.txt'),
        pico_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico_cov.txt'),
        truseq_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq.txt'),
        truseq_cov=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq_cov.txt')
    shell:
        '''
        cd R
        Rscript process_sarantopoulou_create_expression_and_covariates.R {input.v4} {output.v4_mat} {output.v4_cov}
        Rscript process_sarantopoulou_create_expression_and_covariates.R {input.pico} {output.pico_mat} {output.pico_cov}
        Rscript process_sarantopoulou_create_expression_and_covariates.R {input.truseq} {output.truseq_mat} {output.truseq_cov}
        '''
        
rule download_mm9_gtf:
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.gtf.gz')
    shell:
        '''
        wget -O {output} ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz
        '''
        
rule get_mm9_name_id_conversion_table:
    input:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.gtf.gz')
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.id_converstion_table.txt')
    shell:
        '''
        zcat {input} | grep gene_id | perl -pe 's/.+gene_id "(.+?)".+gene_name "(.+?)".+/\\1\\t\\2/g' | sort | uniq > {output}.temp
        cat <(echo -e "gene_id\tgene_name") {output}.temp > {output} 
        rm {output}.temp
        '''
        
rule create_sarantopoulou_matrices_tpm:
    input:
        v4_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4.txt'),
        pico_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico.txt'),
        truseq_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq.txt'),
        gtf=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Mus_musculus.NCBIM37.67.gtf.gz')
    output:
        v4_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4_tpm.txt'),
        pico_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico_tpm.txt'),
        truseq_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq_tpm.txt')
    shell:
        '''
        cd R
        Rscript get_TPM.R {input.v4_mat} {input.gtf} {output.v4_mat}
        Rscript get_TPM.R {input.pico_mat} {input.gtf} {output.pico_mat}
        Rscript get_TPM.R {input.truseq_mat} {input.gtf} {output.truseq_mat}
        '''
        
rule create_sarantopoulou_matrices_tmm:
    input:
        v4_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4.txt'),
        pico_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico.txt'),
        truseq_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq.txt')
    output:
        v4_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4_tmm.txt'),
        pico_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico_tmm.txt'),
        truseq_mat=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq_tmm.txt')
    shell:
        '''
        cd R
        Rscript get_tmm.R {input.v4_mat} {output.v4_mat}
        Rscript get_tmm.R {input.pico_mat} {output.pico_mat}
        Rscript get_tmm.R {input.truseq_mat} {output.truseq_mat}
        '''

rule create_sarantopoulou_lookup_table:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4_tmm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_v4_tpm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico_tpm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_pico_tmm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq_tmm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_truseq_tpm.txt')
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_sarantopoulou.txt')
    shell:
        '''
        cd R
        Rscript process_sarantopoulou_create_lookup_table.R {output} {input}
        '''

##################################################
## Nagarajan

rule download_nagarajan:
    output:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'nagarajan', 'arid1a_knockout_COUNTS.tsv')
    shell:
        '''
        wget -O {output} https://jclab.cruk.cam.ac.uk/share/NG-A51066R2/arid1a_knockout_COUNTS.tsv
        '''
        
rule create_nagaraja_matrix_counts:
    input:
        os.path.join(config['projectDir'], config['expression_datasets_dir'], 'nagarajan', 'arid1a_knockout_COUNTS.tsv')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_cov.txt')
    shell:
        '''
        cd R
        Rscript process_nagarajan_create_expression_and_covariates.R {input} {output}
        '''
        
rule download_grch38_gtf:
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz')
    shell:
        '''
        wget -O {output} ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz
        '''
        
rule create_nagarajan_matrix_tpm:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan.txt'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_tpm.txt')
    shell:
        '''
        cd R
        Rscript get_TPM.R {input} {output}
        '''
        
rule create_nagarajan_matrix_tmm:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan.txt')
    output:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_tmm.txt')
    shell:
        '''
        cd R
        Rscript get_tmm.R {input} {output}
        '''
        
rule create_nagarajan_lookup_table:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_tpm.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_tmm.txt')
    output:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], 'expression_datasets_nagarajan.txt')
    shell:
        '''
        cd R
        Rscript process_nagarajan_create_lookup_table.R {output} {input}
        '''
