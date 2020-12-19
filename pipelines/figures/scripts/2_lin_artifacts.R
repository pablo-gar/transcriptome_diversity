library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
library('purrr')
library('broom')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('../../../R/transcriptome_diversity_tools.R')
source('plot_opts.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Files

transcriptome_diversity_file <- '../data/lin/lin_tpm_transcriptome_diversity.txt'
sample_metadata_file <- '../data/lin/lin_covariates.txt'

#transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations_from_raw/transcriptome_diversity_non_lowly_expressed/lin.txt'
#transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations_from_raw/transcriptome_diversity_non_lowly_expressed/lin_expressed_in_all.txt'

# output
out_prefix <- '../figures/figure_3_features_associated/lin_replicates_'

p_file_strains <- paste0(out_prefix, 'strains.pdf')
p_file_replicates <- paste0(out_prefix, 'replicates.pdf')
p_file_fly_number <- paste0(out_prefix, 'fly_number.pdf')
p_file_sex <- paste0(out_prefix, 'sex.pdf')
p_file_read_depth <- paste0(out_prefix, 'read_depth.pdf')
p_file_model_pvalues <- paste0(out_prefix, 'model_pvalues.pdf')
p_file_combined_final <- paste0(out_prefix, 'final_combined_features.pdf')



#-----------------------------------------------------------------------#
# Read in 

transcriptome_diversity  <- read_tsv(transcriptome_diversity_file)

transcriptome_diversity <- transcriptome_diversity %>% 
    mutate(transcriptome_diversity_rankit=rankitNormalize_vector(transcriptome_diversity))

#Read metadata
sample_metadata <- read_tsv(sample_metadata_file) %>%
    pivot_longer(-ID, names_to='sample_id', values_to='metadata_value') %>%
    pivot_wider(names_from=ID, values_from=metadata_value) %>%
    mutate(mapped_reads=log10(as.double(mapped_reads)))

transcriptome_diversity <- left_join(transcriptome_diversity, sample_metadata)
transcriptome_diversity$dummy <- 'a'
transcriptome_diversity$mapped_reads_rankit <- rankitNormalize_vector(transcriptome_diversity$mapped_reads)


#-----------------------------------------------------------------------#
# Plots

#Plot strain effects
p_strains <- ggplot(transcriptome_diversity, aes(x=DGRP_Number, y=transcriptome_diversity)) +
    geom_violin(fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
    geom_jitter(colour='black', size=p.point_size_tiny, alpha=0.6) +
    theme_grid_y() +
    do.call(theme, c(theme_pars, list(axis.text.x=element_text(angle=30, hjust=1))))

#Plot strain environment effects
label <- transcriptome_diversity %>%
    select(Environment) %>%
    unique() %>%
    mutate(x=0.3, y=Inf, label=paste0('Biological replicate: ', Environment))

minus <- which(names(theme_pars) =='strip.text')
p_environ <- ggplot(transcriptome_diversity, aes(x=transcriptome_diversity)) +
    geom_histogram(fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
    ylab('Count') +
    xlab('Transcriptome diversity') +
    coord_cartesian(xlim=c(min(transcriptome_diversity$transcriptome_diversity), max(transcriptome_diversity$transcriptome_diversity))) +
    facet_grid(Environment~.) +
    geom_text(aes(x=x, y=y, label=label), data=label, fontface='italic', hjust=0, vjust=1, size=p.label_size)  +
    theme_grid_y() +
    do.call(theme, c(theme_pars[-minus], list(strip.text=element_blank())))



#Plot fly number effects
p_fly <- ggplot(transcriptome_diversity, aes(x=transcriptome_diversity)) +
    geom_histogram(fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
    xlab('Transcriptome diversity') +
    ylab('Count') +
    coord_cartesian(xlim=c(min(transcriptome_diversity$transcriptome_diversity), max(transcriptome_diversity$transcriptome_diversity))) +
    facet_grid(Fly_Number~.) +
    theme_grid_y() +
    do.call(theme, c(theme_pars))

#Plot sex effects
label <- transcriptome_diversity %>%
    select(Sex) %>%
    unique() %>%
    mutate(x=0.3, y=Inf, label=paste0('Sex: ', ifelse(Sex=='M', 'Male', 'Female')))

p_sex <- ggplot(transcriptome_diversity, aes(x=transcriptome_diversity)) +
    geom_histogram(fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
    xlab('Transcriptome diversity') +
    ylab('Count') +
    coord_cartesian(xlim=c(min(transcriptome_diversity$transcriptome_diversity), max(transcriptome_diversity$transcriptome_diversity))) +
    geom_text(aes(x=x, y=y, label=label), data=label, fontface='italic', hjust=0, vjust=1, size=p.label_size)  +
    facet_grid(Sex~.) +
    theme_grid_y() +
    do.call(theme, c(theme_pars[-minus], list(strip.text=element_blank())))



#Plot read depth
p_depth_rankit <- scatter(as.data.frame(transcriptome_diversity), x='mapped_reads_rankit', y='transcriptome_diversity_rankit', 
                          method_cor='pearson', regression=T, alpha=0.4, pSize=p.point_size_small, labelSize=p.label_size) +
    xlab('log10(Mapped reads)') +
    ylab('Transcriptome diversity') +
    theme_noGrid() +
    do.call(theme, theme_pars)

p_depth <- scatter(as.data.frame(transcriptome_diversity), x='mapped_reads', y='transcriptome_diversity', 
                   method_cor='spearman',regression=T, alpha=0.4, pSize=p.point_size_small, labelSize=p.label_size) +
    xlab('log10(Mapped reads)') +
    ylab('Transcriptome diversity') +
    theme_noGrid() +
    do.call(theme, theme_pars)



## Modeling transcriptome diversity

model_results <- transcriptome_diversity %>%
    group_by(dummy) %>%
    nest() %>%
    mutate(lm_results=map(data, ~ lm(transcriptome_diversity_rankit ~ Environment + Sex + DGRP_Number + mapped_reads_rankit, data=.x)),
           r_2=map(lm_results, ~ paste0('   r^2 = ', signif(summary(.x)$r.squared), 2)),
           lm_results=map(lm_results, tidy)) %>%
    unnest(c(lm_results, r_2)) %>%
    mutate(term=factor(term, levels=term, ordered=T)) %>%
    ungroup() %>%
    mutate(FDR=p.adjust(p.value, method='BH')) %>%
    mutate(term=gsub('_NumberDGRP', '', term)) %>%
    mutate(term=gsub('mapped_reads_rankit', 'Read depth', term)) %>%
    mutate(term=gsub('SexM', 'Sex(M vs F)', term)) %>%
    mutate(term=gsub('Environment2', 'Replicate(vs 2)', term)) %>%
    mutate(term=gsub('Environment3', 'Replicate(vs 3)', term)) 
    

#Plot
p_model_results <- ggplot(model_results, aes(x=term, y=-log10(FDR))) +
    geom_bar(stat='identity', fill=p.bar_gray_fill, colour=p.bar_gray_colour) +
    geom_text(aes(label=r_2), x=-Inf, y=Inf, data=unique(model_results[,c('r_2')]), vjust=1, hjust=0) +
    geom_hline(yintercept=1.3, linetype='dashed', colour=p.color_dark_blue) +
    annotate(geom='text', x=2, y=1.6, label='FDR = 0.05', fontface='italic', colour=p.color_dark_blue, size=p.label_size, vjust=0, hjust=0) +
    xlab('Regression term') +
    theme_grid_y() +
    do.call(theme, c(theme_pars, list(axis.text.x=element_text(angle=75, hjust=1))))


#Visualize effects on transcriptome diversity of sex, Genotype, and environment
medians <- transcriptome_diversity %>%
    group_by(Environment) %>%
    summarise(yintercept=median(transcriptome_diversity_rankit)) %>%
    ungroup()

label <- transcriptome_diversity %>%
    select(Environment) %>%
    unique() %>%
    mutate(x=-Inf, y=Inf, label=paste0('Biological replicate: ', Environment))

minus <- which(names(theme_pars) =='strip.text')

p_final <- ggplot(transcriptome_diversity, aes(x=log10(mapped_reads), y=transcriptome_diversity_rankit)) +
    geom_hline(aes(yintercept=yintercept), data=medians, linetype='dashed') +
    geom_point(aes(colour=Sex), size=p.point_size_small) +
    geom_text(aes(x=x, y=y, label=label), data=label, fontface='italic', hjust=0, vjust=1, size=p.label_size)  +
    xlab('log10(Mapped reads)') +
    ylab('Transcriptome diversity') +
    facet_wrap(~Environment, nrow=1) +
    scale_fill_manual(values=c('grey40', 'grey60', 'grey80')) +
    scale_colour_manual(values=c(p.negative_color, p.positive_color)) +
    theme_noGrid() +
    do.call(theme, c(theme_pars[-minus], list(strip.text=element_blank())))

ggsave(p_file_strains, p_strains, width=5, height=2)
ggsave(p_file_replicates, p_environ, width=3, height=3)
ggsave(p_file_fly_number, p_fly, width=3, height=5)
ggsave(p_file_sex, p_sex, width=3, height=3)
ggsave(p_file_read_depth, p_depth, width=2.5, height=2.5)
#ggsave, p_depth_rankit, width=1.5, height=1,5)
ggsave(p_file_model_pvalues, p_model_results, width=6, height=3)
ggsave(p_file_combined_final, p_final, width=5, height=2.3)
