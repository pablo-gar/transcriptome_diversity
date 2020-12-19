library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
library('stringr')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('../../../R/misc.R')
source('plot_opts.R') # Contains formating options for all plots
source('gtex_utils.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Input
eqtl_results_dir <- '../data/keele/eqtl/'
eqtl_results_dir <- '/scratch/users/paedugar/transcriptome_diversity/eqtl_results/keele'
eqtl_batch_folder <- file.path(eqtl_results_dir, 'batch')
eqtl_trans_batch_folder <- file.path(eqtl_results_dir, 'transcriptome_diversity+batch')

# output
out_prefix <- '../figures/figure_5_eqtl_results/keele_'

p_file_scatter <- paste0(out_prefix, 'scatter.pdf')

#-----------------------------------------------------------------------#
# Read in 

read_eqtl <- function(eqtl_dir, covariate_type) {
    
    files <- list.files(eqtl_dir, full.names=T, recursive=T, pattern='txt$')
    ids <- str_split_fixed(basename(files), '\\.', n=2)[,1]
    cors <- concatenate_table_files2(files, id_names=ids, delim='\t', header=T, read_function=read_delim, progress=F, col=cols())
    
    cors$nominal_pvalue <- -log10(cors$nominal_pvalue)
    colnames(cors)[colnames(cors) == 'nominal_pvalue'] <- paste0('pvalue_', covariate_type)
    
    return(cors)
}

eqtl_batch <- read_eqtl(eqtl_batch_folder, 'batch_as_covariates') 
eqtl_trans_batch <- read_eqtl(eqtl_trans_batch_folder, 'transcriptome_diversity_batch_as_covariates')

eqtl_results <- left_join(eqtl_batch, eqtl_trans_batch, by=c('gene_id', 'id_names'))

#-----------------------------------------------------------------------#
# plot

p_scatter <- ggplot(eqtl_results, aes(x=pvalue_batch_as_covariates, y=pvalue_transcriptome_diversity_batch_as_covariates)) +
    geom_point(alpha=0.3, size=p.point_size_tiny) +
    #geom_bin2d(bins=100) +
    xlab('eQTL -log10(p-value)\ny = SNP + Batch') +
    ylab('eQTL -log10(p-value)\ny = SNP + Batch + Transcriptome diversity') +
    coord_cartesian(xlim=c(0, 14), ylim=c(0,14)) +
    facet_grid(.~id_names) +
    geom_abline(slope=1, intercept=0, colour=p.abline_colour) +
    theme_noGrid() +
    do.call(theme, theme_pars) 

p_scatter

ggsave(p_file_scatter, p_scatter, height=2.5, width=6)
