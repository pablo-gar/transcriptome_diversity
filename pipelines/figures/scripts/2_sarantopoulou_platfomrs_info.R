library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
library('purrr')
library('broom')
library('stringr')
source('../../../R/plots.R')
source('../../../R/misc.R')
source('../../../R/ggthemes.R')
source('../../../R/transcriptome_diversity_tools.R')
source('plot_opts.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Files
working_dir <- '../data/sarantopoulou/'
out_prefix <- '../figures/2_sarantopoulou_platforms_'

#-----------------------------------------------------------------------#
# Read in 
read_files <- function(x, pattern, ...) {
    
    current_files <- list.files(x, full.names=T, recursive=T, pattern=pattern)
    current_ids <- str_remove(basename(current_files), '.txt')
    current <- concatenate_table_files2(current_files, id_names=current_ids, ...)
    
    return(current)
}


transcriptome_diversity <- read_files(working_dir, 'trans')
transcriptome_diversity$Platform <- str_split_fixed(transcriptome_diversity$id_names, '_', 3)[,2]
transcriptome_diversity$id_names <- NULL

sample_metadata <- read_files(working_dir, 'cov', join_rows=F, join_cols=F) %>%
    map_dfr(function(x) pivot_longer(x, cols=-all_of(c('ID', 'id_names')), names_to='sample_id', values_to='metadata_value')) %>% 
    select(-id_names) %>%
    pivot_wider(sample_id, names_from=ID, values_from=metadata_value)

transcriptome_diversity <- left_join(transcriptome_diversity, sample_metadata)
transcriptome_diversity$mapped_reads <- as.numeric(transcriptome_diversity$mapped_reads)
transcriptome_diversity$log10_mapped_reads <- log10(transcriptome_diversity$mapped_reads)
transcriptome_diversity$mapped_reads_rankit <- rankitNormalize_vector(transcriptome_diversity$mapped_reads)
transcriptome_diversity$transcriptome_diversity_rankit <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
transcriptome_diversity$strain <- factor(transcriptome_diversity$strain, levels=rev(unique(transcriptome_diversity$strain)), ordered=T)


#-----------------------------------------------------------------------#
# Plots

#Plot individual effects
toPlot <- transcriptome_diversity %>%
    mutate(individual = factor(individual, ordered=T, levels=names(sort(tapply(transcriptome_diversity, individual, median)))))

p_individual <- ggplot(toPlot, aes(x=individual, y=transcriptome_diversity)) +
    geom_bar(aes(fill=Platform), colour='black', position='dodge', stat='identity') +
    theme_grid_y() +
    theme(axis.text.x=element_text(angle=30, hjust=1))

#Plot treatment effects
toPlot <- transcriptome_diversity %>%
    mutate(Platform = factor(Platform, ordered=T, levels=names(sort(tapply(transcriptome_diversity, Platform, median)))))

toPlot_means <- toPlot %>%
    group_by(Platform, strain) %>%
    summarise(transcriptome_diversity=mean(transcriptome_diversity)) %>%
    ungroup()

p_plat_unt <- ggplot(filter(toPlot, strain=='UNT'), aes(x=Platform, y=transcriptome_diversity)) +
    geom_violin(fill=p.bar_gray_fill) +
    geom_jitter(data=filter(toPlot, strain=='UNT'), colour='black', width= 0.3, size=p.point_size_medium) +
    ylab('Transcriptome diversity') +
    theme_grid_y() +
    theme(axis.text.x=element_text(angle=30, hjust=1))

p_plat_ilb <- ggplot(filter(toPlot, strain=='ILB'), aes(x=Platform, y=transcriptome_diversity)) +
    geom_violin(fill=p.bar_gray_fill) +
    geom_jitter(data=filter(toPlot, strain=='ILB'), colour='black', width= 0.2, size=p.point_size_medium) +
    ylab('Transcriptome diversity') +
    theme_grid_y() +
    theme(axis.text.x=element_text(angle=30, hjust=1))

#Plot treatment plus Platform effects
p_plat_treatment <- ggplot(transcriptome_diversity, aes(x=strain, y=transcriptome_diversity)) +
    geom_violin(fill=p.bar_gray_fill) +
    geom_jitter(colour='black', width= 0.1, size=0.8, alpha=0.6) +
    facet_grid(.~Platform) +
    ylab('Transcriptome diversity') +
    xlab('Treatment') +
    coord_cartesian(ylim=c(0.62, 0.68)) +
    theme_grid_y() +
    do.call(theme, theme_pars)

#Plot read depth relationship with Platform
p_depth_rankit <- scatter(as.data.frame(transcriptome_diversity), x='mapped_reads_rankit', y='transcriptome_diversity_rankit', 
                          method_cor='pearson', regression=T, alpha=1, pSize=p.point_size_medium, labelSize=p.label_size, facet_x='Platform') +
    xlab('log10(Mapped reads)') +
    ylab('Transcriptome diversity') +
    theme_noGrid() +
    do.call(theme, theme_pars)

p_depth <- scatter(as.data.frame(transcriptome_diversity), x='log10_mapped_reads', y='transcriptome_diversity', 
                   method_cor='spearman',regression=T, alpha=1, pSize=p.point_size_medium, labelSize=p.label_size, facet_x='Platform') +
    xlab('log10(Mapped reads)') +
    ylab('Transcriptome diversity') +
    theme_noGrid() +
    do.call(theme, theme_pars)



ggsave(paste0(out_prefix, 'individual.pdf'), p_individual, width=5, height=3)
ggsave(paste0(out_prefix, 'Platform_untreated.pdf'), p_plat_unt, width=3, height=3)
ggsave(paste0(out_prefix, 'Platform_treat_ILB.pdf'), p_plat_ilb, width=3, height=3)
ggsave(paste0(out_prefix, 'Platform_both_treat_untreat.pdf'), p_plat_treatment, width=3, height=3)
ggsave(paste0(out_prefix, 'read_depth.pdf'), p_depth, width=6, height=2.5)
ggsave(paste0(out_prefix, 'read_depth_rankit.pdf'), p_depth_rankit, width=6, height=2.5)
