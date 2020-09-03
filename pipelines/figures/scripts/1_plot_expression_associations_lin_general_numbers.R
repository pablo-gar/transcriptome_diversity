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
plot_file_N <- '../figures/1_all_expression_associations_lin_N.pdf'
plot_file_directionality_tpm <- '../figures/1_all_expression_associations_lin_directionality_tpm.pdf'
plot_file_directionality_tmm <- '../figures/1_all_expression_associations_lin_directionality_tmm.pdf'

# Files
association_results_tpm <- '../data/lin_expression_associations_tpm.txt'
association_results_tmm <- '../data/lin_expression_associations_tmm.txt'


FDR_cutoff <- 0.05 


#-----------------------------------------------------------------------#
# Read in 

expression_associations_tpm <- read_tsv(association_results_tpm)
expression_associations_tmm <- read_tsv(association_results_tmm)



#-----------------------------------------------------------------------#
# Process data 

expression_associations_tpm$Type <- 'TPM'
expression_associations_tmm$Type <- 'TMM'

# Merges tables
expression_associations <- bind_rows(expression_associations_tmm, expression_associations_tpm) %>%
    mutate(Type=factor(Type, levels=c('TPM', 'TMM'), ordered=T))
                                     
#-----------------------------------------------------------------------#
# Calculations

significant_associations <- expression_associations %>%
    group_by(Type) %>%
    summarise(total =n(), significant=sum(FDR < FDR_cutoff)) %>%
    ungroup() %>%
    mutate(non_significant= total - significant, Significant = 100*significant/total, Non_significant=100-Significant,
            ) 

# Percentages for barplot
significant_associations_for_bar <- significant_associations %>%
    pivot_longer(cols=c(Significant, Non_significant), names_to='FDR', values_to='Percentage')

# Absolute numbers for lables on plot
lab <- significant_associations %>%
    pivot_longer(cols=c(significant, non_significant), names_to='FDR', values_to='N')  %>%
    mutate(lab=paste0('n = ', N))
lab$Percentage <- lab$Significant
lab$Percentage[lab$FDR == 'significant'] <- 0


# number of positive and negative associations
associations_sign <- expression_associations %>%
    filter(FDR < FDR_cutoff) 

lab_associations <- associations_sign %>%
    group_by(Type) %>%
    summarise(lab=paste0('\n      positive = ', sum(statistic>0), '\n',
                         '      negative = ', sum(statistic<0))) %>%
    ungroup() %>%
    mutate(x=-Inf, y=Inf)
    

#-----------------------------------------------------------------------#
# plots

# Number of significant and non-significant associations
p_n_associations <- ggplot(significant_associations_for_bar, aes(y=Type, x=Percentage)) +
    geom_bar(aes(fill=FDR), width=0.8, colour='black', stat='identity') +
    geom_text(aes(label=lab, x=Percentage + 2), colour='white', data=lab, hjust=0, fontface='bold.italic', size=p.label_size-0.3) +
    scale_fill_manual(values=c(p.negative_color, p.positive_color)) +
    theme_noGrid() +
    do.call(theme, c(theme_pars, list(legend.position='top'))) +
    do.call(guides, guide_pars)


# Directionality of significant associations
p_directionality_tpm <- ggplot(filter(associations_sign, Type == 'TPM'), aes(statistic)) +
    geom_histogram(fill=p.bar_gray_fill, colour=p.bar_gray_colour, bins=40) +
    geom_text(aes(label=lab, x=x, y=y), data=filter(lab_associations, Type=='TPM'), hjust=0, vjust=1, fontface='italic', size=p.label_size) +
    xlab('F-statistic\n(Directionality of association)') + 
    ylab('Count') +
    theme_grid_x() +
    do.call(theme, theme_pars) 

p_directionality_tmm <- ggplot(filter(associations_sign, Type == 'TMM'), aes(statistic)) +
    geom_histogram(fill=p.bar_gray_fill, colour=p.bar_gray_colour, bins=40) +
    geom_text(aes(label=lab, x=x, y=y), data=filter(lab_associations, Type=='TMM'), hjust=0, vjust=1, fontface='italic', size=p.label_size) +
    xlab('F-statistic\n(Directionality of association)') + 
    ylab('Count') +
    theme_grid_x() +
    do.call(theme, theme_pars) 

ggsave(plot_file_N, p_n_associations, height=1.5, width=4)
ggsave(plot_file_directionality_tpm, p_directionality_tpm, height=1.610, width=4)
ggsave(plot_file_directionality_tmm, p_directionality_tmm, height=3, width=3)

