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

FDR_cutoff <- 0.05 

# Input
expression_associations_outdir <- '../data/gtex/expression_associations/tpm'
expression_associations_outdir_tmm <- '../data/gtex/expression_associations/tmm'

PCA_out_dir_gtex <- '../data/gtex/PCA_var_explained_by_transcriptome_diversity/tpm'
PCA_out_dir_gtex_tmm <- '../data/gtex/PCA_var_explained_by_transcriptome_diversity/tmm'

# output
out_prefix <- '../figures/figure_4_big_consortia/gtex_expression_associations_'

p_file_tpm <- paste0(out_prefix, 'tpm.pdf')
p_file_tmm <- paste0(out_prefix, 'tmm.pdf')

out_prefix_PCs <- '../figures/figure_4_big_consortia/gtex_PC_cors_'

p_file_tpm_PCs <- paste0(out_prefix_PCs, 'tpm.pdf')
p_file_tmm_PCs <- paste0(out_prefix_PCs, 'tmm.pdf')

#-----------------------------------------------------------------------#
# PCA 
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#
# Read in 
read_files <- function(x, pattern, type, ...) {
    
    current_files <- list.files(x, full.names=T, recursive=T, pattern=pattern)
    current_ids <- str_remove(basename(dirname(current_files)), '.txt')
    current_ids <- str_remove(basename(current_ids), 'gtex_')
    current <- concatenate_table_files2(current_files, id_names=current_ids, col=cols(), ...)
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

# Clean up tissues
pca_all <- filter(pca_all, !id_names %in% GTEX_TISSUES_LOW_SAMPLES)
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

p_tpm_PCs <- p_cor_PC1_tpm + p_signs_tpm + p_top_PC_tpm + p_top_PC_var_explained_tpm + plot_layout(design="AAAAAAAAABCD")
p_tmm_PCs <- p_cor_PC1_tmm + p_signs_tmm + p_top_PC_tmm + p_top_PC_var_explained_tmm + plot_layout(design="AAAAAAAAABCD")

ggsave(p_file_tpm, p_tpm_PCs, width=4.8, height=7.2)
ggsave(p_file_tmm, p_tmm_PCs, width=4.8, height=7.2)


#-----------------------------------------------------------------------#
# EXPRESSION ASSOCIATIONS
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#
# Read in 
read_files <- function(x, pattern, type, ...) {
    
    current_files <- list.files(x, full.names=T, recursive=T, pattern=pattern)
    current_ids <- str_remove(basename(current_files), '.txt')
    current_ids <- str_remove(basename(current_ids), 'gtex_')
    current <- concatenate_table_files2(current_files, id_names=current_ids, ...)
    current$Type <- type
    
    return(current)
}

exp_associations_tpm <- read_files(expression_associations_outdir, '.txt', 'TPM', col_types = cols())
exp_associations_tmm <- read_files(expression_associations_outdir_tmm, '.txt', 'TMM', col_types = cols())

#Merge
expression_associations <- bind_rows(exp_associations_tpm, exp_associations_tmm)

# Remove tissues with low number of samples
expression_associations <- filter(expression_associations, !id_names %in% GTEX_TISSUES_LOW_SAMPLES)
# Correct names for readability
expression_associations$id_names <- correct_tissue_names(expression_associations$id_names)



#-----------------------------------------------------------------------#
# Process

# No. significant associations
significant_associations <- expression_associations %>%
    group_by(Type, id_names) %>%
    summarise(total =n(), significant=sum(FDR < FDR_cutoff)) %>%
    ungroup() %>%
    mutate(non_significant= total - significant, Significant = 100*significant/total, Non_significant=100-Significant) 

# Percentages for barplot
significant_associations_for_bar <- significant_associations %>%
    pivot_longer(cols=c(Significant, Non_significant), names_to='FDR', values_to='Percentage') 

# Positive vs negative
pos_neg <- expression_associations %>%
    filter(FDR < FDR_cutoff) %>% 
    group_by(Type, id_names) %>%
    summarise(pos=sum(statistic>0), neg=sum(statistic<0)) %>%
    ungroup() %>%
    mutate(dummy=1)


#-----------------------------------------------------------------------#
# plots

# Number of significant and non-significant associations
plot_n_associations <- function(x) {
    
    ggplot(x, aes(y=id_names, x=Percentage)) +
        geom_bar(aes(fill=FDR), width=0.8, colour='black', stat='identity') +
        scale_fill_manual(values=c(p.negative_color, p.positive_color)) +
        xlab("Percentage of \n genes tested") +
        ylab("") +
        geom_vline(xintercept=25, linetype='dotted') +
        geom_vline(xintercept=50, linetype='dotted') +
        geom_vline(xintercept=75, linetype='dotted') +
        theme_noGrid() +
        do.call(theme, c(theme_pars, list(legend.position='top', axis.text.y=element_blank()))) +
        do.call(guides, guide_pars)
    
}

# To plot number of positive vs negative associations
plot_text_column <- function(x, label, colour) {
    ggplot(x, aes(y=id_names, x=dummy)) +
        geom_text(aes_string(label=label), colour=colour, size=p.label_size) +
        #theme_bw()
        theme_fullBlank()
}


# Separate between expression types 
signif_tpm <- filter(significant_associations_for_bar, Type=='TPM')
signif_tmm <- filter(significant_associations_for_bar, Type=='TMM')
pos_neg_tpm <- filter(pos_neg, Type=='TPM')
pos_neg_tmm <- filter(pos_neg, Type=='TMM')

# Order tissues based on % of significant associations
#signif_tpm <- mutate(signif_tpm, id_names=factor(id_names, levels=unique(id_names[order(significant/total)])))
#signif_tmm <- mutate(signif_tmm, id_names=factor(id_names, levels=unique(id_names[order(significant/total)])))
signif_tpm <- mutate(signif_tpm, id_names=factor(id_names, levels=levels(top_PC_tpm$id_names)))
signif_tmm <- mutate(signif_tmm, id_names=factor(id_names, levels=levels(top_PC_tmm$id_names)))
pos_neg_tpm <- mutate(pos_neg_tpm, id_names=factor(id_names, levels=levels(signif_tpm$id_names)))
pos_neg_tmm <- mutate(pos_neg_tmm, id_names=factor(id_names, levels=levels(signif_tmm$id_names)))

# plots
p_n_associations_tpm <- plot_n_associations(x=signif_tpm)
p_n_associations_tmm <- plot_n_associations(x=signif_tmm)

p_pos_tpm <- plot_text_column(x=pos_neg_tpm, label='pos', colour=p.positive_color_green)
p_pos_tmm <- plot_text_column(x=pos_neg_tmm, label='pos', colour=p.positive_color_green)

p_neg_tpm <- plot_text_column(x=pos_neg_tpm, label='neg', colour=p.negative_color_red)
p_neg_tmm <- plot_text_column(x=pos_neg_tmm, label='neg', colour=p.negative_color_red)

# final plots
p_tpm <- p_n_associations_tpm + p_pos_tpm + p_neg_tpm + plot_layout(design="AAAAAAAAABC")
p_tmm <- p_n_associations_tmm + p_pos_tmm + p_neg_tmm + plot_layout(design="AAAAAAAAABC")


ggsave(p_file_tpm_PCs, p_tpm, width=2.65, height=7.2)
ggsave(p_file_tmm_PCs, p_tmm, width=2.65, height=7.2)

