library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('plot_opts.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Outputs
plot_file_var_explained_prefix <- '../figures/figure_2_example/lin_pcs_'

plot_file_PCs_tpm <- paste0(plot_file_var_explained_prefix, 'PC1_PC2_tpm.pdf')
plot_file_PCs_tmm <- paste0(plot_file_var_explained_prefix, 'PC1_PC2_tmm.pdf')
plot_file_var_explained_tpm <- paste0(plot_file_var_explained_prefix, 'var_explained_tpm.pdf')
plot_file_var_explained_tmm <- paste0(plot_file_var_explained_prefix, 'var_explained_tmm.pdf')
plot_file_cors_trans_tpm <- paste0(plot_file_var_explained_prefix, 'cor_trans_tpm.pdf')
plot_file_cors_trans_tmm <- paste0(plot_file_var_explained_prefix, 'cor_trans_tmm.pdf')


# Files
covs_file <- '../data/lin/lin_covariates.txt'

cors_tmm_file <- '../data/lin/lin_tmm_pca/cors.txt'
stats_tmm_file <- '../data/lin/lin_tmm_pca/PC_stats.txt'
pcs_tmm_file <- '../data/lin/lin_tmm_pca/PC_matrix.txt'

cors_tpm_file <- '../data/lin/lin_tpm_pca/cors.txt'
stats_tpm_file <- '../data/lin/lin_tpm_pca/PC_stats.txt'
pcs_tpm_file <- '../data/lin/lin_tpm_pca/PC_matrix.txt'

# PCs to plot for var explained plots
pcs_to_plot <- 8


#-----------------------------------------------------------------------#
# Read in 

covs <- read_tsv(covs_file) %>% pivot_longer(-ID, names_to='sample_id', values_to='Sex') %>% filter(ID=='Sex')

cors_tmm <- read_tsv(cors_tmm_file) %>% mutate(Type='TMM')
stats_tmm <- read_tsv(stats_tmm_file) %>% mutate(Type='TMM') %>% pivot_longer(-c(stat,Type), names_to='feature', values_to='value')
pcs_tmm <- read_tsv(pcs_tmm_file) %>% mutate(Type='TMM')

cors_tpm <- read_tsv(cors_tpm_file) %>% mutate(Type='TPM') 
stats_tpm <- read_tsv(stats_tpm_file) %>% mutate(Type='TPM') %>% pivot_longer(-c(stat,Type), names_to='feature', values_to='value')
pcs_tpm <- read_tsv(pcs_tpm_file) %>% mutate(Type='TPM')


stats <- bind_rows(stats_tpm, stats_tmm) 
cors <- bind_rows(cors_tpm, cors_tmm) 
pcs <- bind_rows(select(pcs_tmm, c(PC1, PC2, sample_id, transcriptome_diversity, Type)), select(pcs_tpm, c(PC1, PC2, sample_id, transcriptome_diversity, Type))) %>%
    mutate(Type=factor(Type, levels=c('TPM', 'TMM'), ordered=T)) %>% 
    left_join(covs)

#-----------------------------------------------------------------------#
# Process data 

pc_var_explained_cor_trans <- stats %>%
    filter(stat=='Proportion of Variance') %>%
    mutate(value=value*100) %>% 
    left_join(cors) %>%
    mutate(PC=as.numeric(gsub('PC', '', feature)), pvalue_bonf=-log10(pvalue_bonf)) %>% 
    filter(PC <= pcs_to_plot) %>%
    dplyr::rename(var_explained='value') %>%
    select(Type, PC, var_explained, pvalue_bonf, pearson_cor) %>% 
    mutate(pearson_cor=abs(pearson_cor))
    

#-----------------------------------------------------------------------#
# plot


# PCs
scatter <- function(x, type) {
    
    ggplot(filter(x, Type==type), aes(x=PC1, y=PC2)) +
       geom_point(aes(colour=transcriptome_diversity, shape=Sex), size=p.point_size_small) +
       xlab('PC1') +
       ylab('PC2') +
       scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
       #scale_colour_gradient(low = "#307cff", high = "#ffe645", na.value = NA) +
       theme_noGrid() +
       do.call(theme, c(theme_pars, list(legend.position='top')))
}

p_tpm <- scatter(pcs, 'TPM')
p_tmm <- scatter(pcs, 'TMM')


# Var explained and correlations with PCs and transcriptome_diversity
bar <- function(toPlot, type, y, ylabel) {
    
    ggplot(filter(toPlot, Type==type), aes_string(x='PC',y=y)) +
        geom_bar(stat='identity', width=0.7, fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
        xlab('Principal Component') +
        ylab(ylabel) +
        theme_grid_y() +
        do.call(theme, c(theme_pars, list(legend.position='top')))
    
}

p_cors_trans_tpm <- bar(pc_var_explained_cor_trans, 'TPM', 'pearson_cor', ylabel="Transcriptome\n diversity association\n|r|")
p_var_exp_tpm <- bar(pc_var_explained_cor_trans, 'TPM', 'var_explained', ylabel='Variance\n explained (%)')

p_cors_trans_tmm <- bar(pc_var_explained_cor_trans, 'TMM', 'pearson_cor', ylabel="Transcriptome\n diversity association\n|r|")
p_var_exp_tmm <- bar(pc_var_explained_cor_trans, 'TMM', 'pearson_cor', ylabel='Variance\n explained (%)')

ggsave(plot_file_PCs_tpm, p_tpm, height=3.2, width=3)
ggsave(plot_file_PCs_tmm, p_tmm, height=3.2, width=3)
ggsave(plot_file_var_explained_tpm, p_var_exp_tpm, height=1.5, width=3.5)
ggsave(plot_file_var_explained_tmm, p_var_exp_tmm, height=1.5, width=3.5)
ggsave(plot_file_cors_trans_tpm, p_cors_trans_tpm, height=1.5, width=3.5)
ggsave(plot_file_cors_trans_tmm, p_cors_trans_tmm, height=1.5, width=3.5)
