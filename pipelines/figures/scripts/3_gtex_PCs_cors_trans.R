library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
library('stringr')
library('patchwork')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('../../../R/misc.R')
source('../../../R/transcriptome_diversity_tools.R')
source('plot_opts.R') # Contains formating options for all plots
source('gtex_utils.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Input
PCA_out_dir_gtex <- '../data/gtex/PCA_var_explained_by_transcriptome_diversity/tpm'
PCA_out_dir_gtex_tmm <- '../data/gtex/PCA_var_explained_by_transcriptome_diversity/tmm'

# output
out_prefix <- '../figures/figure_4_big_consortia/gtex_PC_cors_'

p_file_tpm <- paste0(out_prefix, 'tpm.pdf')
p_file_tmm <- paste0(out_prefix, 'tmm.pdf')

#-----------------------------------------------------------------------#
# Read in 
read_files <- function(x, pattern, type, ...) {
    
    current_files <- list.files(x, full.names=T, recursive=T, pattern=pattern)
    current_ids <- str_remove(basename(dirname(current_files)), '.txt')
    current <- concatenate_table_files2(current_files, id_names=current_ids, col=cols() ...)
    current$type <- type
    
    return(current)
}

# GTEx
cors_gtex <- read_files(PCA_out_dir_gtex, 'cors.txt', 'TPM')
pca_stats_gtex <-read_files(PCA_out_dir_gtex, 'PC_stats', 'TPM')


# GTEx tmm 
cors_gtex_tmm <- read_files(PCA_out_dir_gtex_tmm,  'cors.txt', 'TMM')
pca_stats_gtex_tmm <-read_files(PCA_out_dir_gtex_tmm,  'PC_stats', 'TMM')

#Merge
cors_all <- bind_rows(cors_gtex, cors_gtex_tmm)
pca_stats_all <- bind_rows(pivot_longer(pca_stats_gtex, starts_with('PC'), names_to='feature', values_to='stat_value'),
                            pivot_longer(pca_stats_gtex_tmm, starts_with('PC'), names_to='feature', values_to='stat_value'))

pca_stats_all <- filter(pca_stats_all, stat=='Proportion of Variance') %>%
    dplyr::rename(proportion_var = 'stat_value') %>%
    select(-stat)

pca_all  <- left_join(cors_all,pca_stats_all)
pca_all$id_names <- correct_tissue_names(pca_all$id_names)

#-----------------------------------------------------------------------#
# Process

# Getting the top PC that correlates with transcriptome diversity
top_PC <- pca_all %>%
    group_by(id_names, type) %>%
    arrange(desc(abs(pearson_cor))) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    group_by(type) %>%
    mutate(dummy='a', cor_sign = ifelse(pearson_cor>0, '+', '-'), proportion_var=signif(proportion_var,3))

top_PC_tpm <- filter(top_PC, type=='TPM')
top_PC_tmm <- filter(top_PC, type=='TMM')

top_PC_tpm$id_names <- factor(top_PC_tpm$id_names, levels=top_PC_tpm$id_names[order(abs(top_PC_tpm$pearson_cor))], ordered=T)
top_PC_tmm$id_names <- factor(top_PC_tmm$id_names, levels=top_PC_tmm$id_names[order(abs(top_PC_tmm$pearson_cor))], ordered=T)


#-----------------------------------------------------------------------#
# plot

plot_cor <- function(x){
    
    ggplot(x, aes(y=id_names, x=abs(pearson_cor))) +
        geom_point() +
        xlim(c(0,1)) +
        xlab("|Person's r|") +
        ylab("") +
        theme_bw() +
        do.call(theme, theme_pars)
    
}

plot_signs <- function(x) {
    ggplot(x, aes(y=id_names, x=dummy)) +
        geom_text(aes(label=cor_sign)) +
        theme_fullBlank()
}

plot_top_PC <- function(x) {
    ggplot(x, aes(y=id_names, x=dummy)) +
        geom_text(aes(label=feature), size=2.1) +
        theme_fullBlank()

}

plot_var_explained <- function(x) {
    
    ggplot(x, aes(y=id_names, x=dummy)) +
        geom_tile(aes(fill=proportion_var)) +
        geom_text(aes(label=proportion_var), size=2.1) +
        scale_fill_gradientn(colors=c('white', 'red')) +
        theme_fullBlank() +
        theme(legend.position='top') +
        do.call(guides, guide_pars)
    
}

p_cor_PC1_tpm <- plot_cor(top_PC_tpm)
p_cor_PC1_tmm <- plot_cor(top_PC_tmm)

p_signs_tpm <- plot_signs(top_PC_tpm)
p_signs_tmm <- plot_signs(top_PC_tmm)

p_top_PC_tpm <- plot_top_PC(top_PC_tpm)
p_top_PC_tmm <- plot_top_PC(top_PC_tmm)

p_top_PC_var_explained_tpm <- plot_var_explained(top_PC_tpm)
p_top_PC_var_explained_tmm <- plot_var_explained(top_PC_tmm)

p_tpm <- p_cor_PC1_tpm + p_signs_tpm + p_top_PC_tpm + p_top_PC_var_explained_tpm + plot_layout(design="AAAAAAAAABCD")
p_tmm <- p_cor_PC1_tmm + p_signs_tmm + p_top_PC_tmm + p_top_PC_var_explained_tmm + plot_layout(design="AAAAAAAAABCD")

ggsave(p_file_tpm, p_tpm, width=4.8, height=7.2)
ggsave(p_file_tmm, p_tmm, width=4.8, height=7.2)
