#Finds PEER covariates associated with transcriptome diversity

library('readr')
library('dplyr')
library('tidyr')
source('../../R/gtex.R', chdir=T)
source('../../R/misc.R', chdir=T)

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    tissue <- cmdArgs[2]
    cor_pvalue_cutoff <- as.numeric(cmdArgs[3])
    gtex_version <- cmdArgs[4]
    out_cors <- cmdArgs[5]
    out_stats <- cmdArgs[6]
    
    if(!gtex_version %in% c('v7', 'v8'))
        stop('Gtex version has to be v7 or v8 ')
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/transcriptome_diversity/Whole_Blood.txt'
    #tissue <- 'Whole_Blood'
    #cor_pvalue_cutoff <- 0.05
    #gtex_version <- 'v8'
    #out_dir <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/peer_associated_with_transcriptome_diversity/Whole_Blood/'
    
    # Read trasncriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    transcriptome_diversity$transcriptome_diversity <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
    transcriptome_diversity$ind <- gtexLongToShort(transcriptome_diversity$gtexId)
    
    # Read PEER
    if(gtex_version=='v8') {
        peer <- readPEERCovariates(tissue, file.path(GTEX_CON$root, GTEX_CON$gtexPeerV8Dir))
    } else {
        peer <- readPEERCovariates(tissue)
    }
    
    peer <- as.data.frame(rankitNormalize(as.matrix(peer), 2))
    peer <- select(peer, starts_with('InferredCov'))
    peer$ind <- rownames(peer)
    peer <- pivot_longer(peer, -ind, names_to='feature', values_to='feature_value')
    
    # Merge
    merged <- inner_join(transcriptome_diversity, peer)

    cors <- merged %>%
        group_by(feature) %>%
        summarise(pearson_cor = cor_test(transcriptome_diversity, feature_value, val='estimate'),
                  pvalue = cor_test(transcriptome_diversity, feature_value, val='p.value'),
                  signed_pvalue = sign(pearson_cor) * -log10(pvalue)) %>%
        ungroup() %>%
        filter(!is.na(pearson_cor)) %>%
        mutate(pvalue_bonf = p.adjust(pvalue)) %>%
        filter(pvalue_bonf < cor_pvalue_cutoff)
    
    stats <- data.frame(total_peer=length(unique(peer$feature)), cor_with_transcriptome_diversity_peer=nrow(cors), cor_pvalue_cutoff=cor_pvalue_cutoff)
    
    # Save results
    write_tsv(cors, out_cors)
    write_tsv(stats, out_stats)
    
    
}

main()
