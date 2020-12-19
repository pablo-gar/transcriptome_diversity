library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('plot_opts.R') # Contains formating options for all plots
source('gtex_utils.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Input
all_results_file <- '../data/gtex/eqtl/eqtl_results.txt'

# output
out_prefix <- '../figures/figure_5_eqtl_results/gtex_'

p_file_histograms <- paste0(out_prefix, 'p_value_histograms.pdf')

#-----------------------------------------------------------------------#
# Read in 

all_results <- read_tsv(all_results_file)

#-----------------------------------------------------------------------#
# Process

all_results <- dplyr::rename(all_results, covariates='group')

# Human-readable categories
#covariate_types <- setNames(c('None', 'metadata', 'PEER', 'PEER + metadata', 'Transcriptome diversity', 'Transcriptome diversity + metadata'),
#                            c('no_covariates', 'no_peer', 'peer', 'peer_plus' , 'transcriptome_diversity_as_covariates', 'transcriptome_diversity_as_covariates_plus_nonPeer'))

covariate_types <- setNames(c('y = SNP', 'y = SNP + Metadata', 'y = SNP + PEER', 'y = SNP + Metadata + PEER', 'y = SNP + Transcriptome diversity', 'y = SNP + Metadata + Transcriptome diversity'),
                            c('no_covariates', 'no_peer', 'peer', 'peer_plus' , 'transcriptome_diversity_as_covariates', 'transcriptome_diversity_as_covariates_plus_nonPeer'))

all_results$covariates <- factor(covariate_types[all_results$covariates], levels=covariate_types, ordered=T)


#-----------------------------------------------------------------------#
# plot

p_histograms <- ggplot(all_results, aes(x=permutation_pvalue)) +
    geom_histogram(bins=50, colour=p.bar_gray_colour, fill=p.bar_gray_fill) +
    xlab('eQTL p-value') +
    ylab('Count') +
    facet_wrap(.~covariates, ncol=2) +
    theme_grid_y() +
    do.call(theme, theme_pars) 

p_histograms

ggsave(p_file_histograms, p_histograms, width=5, height=6)
