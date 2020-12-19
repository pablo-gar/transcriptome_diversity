library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('plot_opts.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

out_folder <- '../figures/figure_2_example/'

# Outputs
plot_file_tpm <- file.path(out_folder, 'association_example_tpm.pdf')
plot_file_tmm <- file.path(out_folder, 'association_example_tmm.pdf')

# Files
exp_mat_file_tpm <- '../data/lin/lin_tpm.txt'
exp_mat_file_tmm <- '../data/lin/lin_tmm.txt'
trans_diversity_file_tpm  <- '../data/lin/lin_tpm_transcriptome_diversity.txt'


# Choose gene of interest
gene <- 'FBgn0013303' # Neurocalcin: Exhibits calcium ion binding activity. Involved in negative regulation of neuronal signal transduction and positive regulation of circadian sleep/wake cycle, sleep
gene_alias <- 'Nca'


#-----------------------------------------------------------------------#
# Read in 

exp_mat_tpm <- read_tsv(exp_mat_file_tpm)
exp_mat_tmm <- read_tsv(exp_mat_file_tmm)
trans_diversity_tpm  <- read_tsv(trans_diversity_file_tpm)



#-----------------------------------------------------------------------#
# Process data 

exp_mat_tpm$type <- 'TPM'
exp_mat_tmm$type <- 'TMM'

# Gets gene of interest and merges tables
exp_mat <- bind_rows(exp_mat_tmm[exp_mat_tmm$gene_id == gene,],
                     exp_mat_tpm[exp_mat_tmm$gene_id == gene,]) %>% 
    pivot_longer(-c(gene_id, type), names_to='sample_id', values_to='expression')

trans_diversity <- trans_diversity_tpm

all_data <- left_join(exp_mat, trans_diversity) %>%
    mutate(type=factor(type, levels=c('TPM', 'TMM'), ordered=T))

#-----------------------------------------------------------------------#
# plot

p_tpm <- scatter(as.data.frame(filter(all_data, type == 'TPM')), x='expression', y='transcriptome_diversity', regression=F, labelSize=p.label_size, pSize=p.point_size_tiny) + 
    xlab(paste0(gene_alias, ' expression (TPM)')) +
    ylab('Transcriptome diversity') +
    coord_cartesian(ylim=c(0.2,0.85)) +
    theme_noGrid() +
    do.call(theme, theme_pars)

p_tmm <- scatter(as.data.frame(filter(all_data, type == 'TMM')), x='expression', y='transcriptome_diversity', regression=F, labelSize=p.label_size, pSize=p.point_size_tiny) + 
    xlab(paste0(gene_alias, ' expression (TMM)')) +
    ylab('Transcriptome diversity') +
    coord_cartesian(ylim=c(0.2,0.85)) +
    theme_noGrid() +
    do.call(theme, theme_pars)

ggsave(plot_file_tpm, p_tpm, height=3, width=3)
ggsave(plot_file_tmm, p_tmm, height=3, width=3)
