library("readr")
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    n_reads <- as.numeric(cmdArgs[2])
    out_exp_mat <- cmdArgs[3]
    
    # Read counts
    exp_mat <- read_tsv(exp_mat_file)
    
    # Sample
    exp_mat <- resize_expression(exp_mat, n=n_reads, method='montecarlo')
    
    #Read gene annotation for gene lengths
    write_tsv(exp_mat, out_exp_mat)
    
    
}

main()
