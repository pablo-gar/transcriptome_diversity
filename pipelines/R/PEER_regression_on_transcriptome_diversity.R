# Calculates how much variance of transcriptome diversity can PEER covariates explain

library('readr')
library('dplyr')
library('tidyr')
library('broom')
source('../../R/gtex.R', chdir=T)
source('../../R/misc.R', chdir=T)

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    tissue <- cmdArgs[2]
    gtex_version <- cmdArgs[3]
    out_file <- cmdArgs[4]
    
    if(!gtex_version %in% c('v7', 'v8'))
        stop('Gtex version has to be v7 or v8 ')
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/transcriptome_diversity/Whole_Blood.txt'
    #tissue <- 'Whole_Blood'
    #gtex_version <- 'v8'
    
    # Read trasncriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    transcriptome_diversity <- transcriptome_diversity %>%
        mutate(transcriptome_diversity=rankitNormalize_vector(transcriptome_diversity),
               ind=gtexLongToShort(gtexId)) %>%
        select(-gtexId)
    
    # Read PEER
    if(gtex_version=='v8') {
        peer <- readPEERCovariates(tissue, file.path(GTEX_CON$root, GTEX_CON$gtexPeerV8Dir))
    } else {
        peer <- readPEERCovariates(tissue)
    }
    
    peer <- as.data.frame(rankitNormalize(as.matrix(peer), 2))
    peer <- select(peer, starts_with('InferredCov'))
    peer$ind <- rownames(peer)
    
    # Merge
    merged <- inner_join(transcriptome_diversity, peer)
    
    # Do linear regressions
    lm_result <- lm(transcriptome_diversity ~ . , data=select(merged, -ind))
    lm_result_tidy <- tidy(lm_result)
    lm_result_tidy$model_r_squared <- summary(lm_result)[['r.squared']]
    
    # Write results
    write_tsv(lm_result_tidy, out_file)

    
}

main()
