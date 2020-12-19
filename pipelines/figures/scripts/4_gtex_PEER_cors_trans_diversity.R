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
cor_file <- '../data/gtex/PEER/Whole_Blood/cors.txt'

# output
out_prefix <- '../figures/figure_4.1_GTEx_PEER/gtex_PEER_'

p_file_peer_cors_bar <- paste0(out_prefix, 'cors.pdf')

fdr_cutoff <- 0.05

#-----------------------------------------------------------------------#
# Read in 

cors <- read_tsv(cor_file)

#-----------------------------------------------------------------------#
# Process

cors$FDR <- ifelse(cors$fdr < fdr_cutoff, '< 0.05', '>= 0.05')
cors$feature <- gsub('InferredCov', '', cors$feature)
cors$feature <- factor(cors$feature, levels=cors$feature[order(cors$fdr)])

#-----------------------------------------------------------------------#
# plot

p_peer_cors_bar <- ggplot(cors, aes(x=feature, y=abs(pearson_cor), fill=FDR)) +
    geom_bar(stat='identity', position='dodge', colour=p.bar_gray_colour) +
    ylab("Correlation with\ntranscriptome diversity\n(Spearman's r)") +
    xlab('PEER covariate') +
    scale_fill_manual(values=c(p.positive_color, p.negative_color)) +
    theme_grid_y() +
    theme(legend.position='top') +
    do.call(theme, theme_pars) 


ggsave(p_file_peer_cors_bar, p_peer_cors_bar, width=6.7, height=2.5)
