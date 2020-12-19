library('ggplot2')
library('readr')
library('dplyr')
library('tidyr')
library('stringr')
source('../../../R/plots.R')
source('../../../R/ggthemes.R')
source('plot_opts.R') # Contains formating options for all plots
source('gtex_utils.R') # Contains formating options for all plots

#-----------------------------------------------------------------------#
# Args

# Inputs
in_dir <- '../data/nagarajan/differential_gene_expression'
in_files <- file.path(in_dir, c('edgeR_results_fulvestrant_wt.Rds', 'edgeR_results_genotype.Rds', 'edgeR_results_jq1_wt.Rds', 'edgeR_results_tamoxifen_wt.Rds'))

# output
out_prefix <- '../figures/figure_6_differential_expression/nagarajan_'

p_file_arid <- paste0(out_prefix, 'wt_vs_arid.pdf')
p_file_fulvestrant <- paste0(out_prefix, 'wt_vs_fulvestrant.pdf')
p_file_jq1 <- paste0(out_prefix, 'wt_vs_jq1.pdf')
p_file_tamoxifen <- paste0(out_prefix, 'wt_vs_tamoxifen.pdf')


#-----------------------------------------------------------------------#
# Read in 

differential <- list()
for(i in in_files){
    
    c_name <- basename(i) %>%
        str_remove( '.Rds') %>%
        str_remove('edgeR_results_')
    
    differential[[c_name]] <- readRDS(i)
}


#-----------------------------------------------------------------------#
# Process

make_differential_table <- function(x) {
    
    final_table <- as_tibble(x[[1]]$differential_genes$table)
    for (i in 2:length(x)) {
       
        c_table <- tibble(a=x[[i]]$differential_genes$table$PValue)
        colnames(c_table) <- paste0('PValue_', names(x)[i])
        final_table <- bind_cols(final_table, c_table)
        
    }

    final_table <- final_table %>%
        pivot_longer(starts_with('PValue_'), names_to='type', values_to='PValue_cov') %>%
        mutate(PValue=-log10(PValue), PValue_cov=-log10(PValue_cov), is_1.2_order_magnitude=ifelse(PValue_cov-PValue > 1.2, 'bigger', ifelse(PValue_cov-PValue < -1.2, 'smaller', 'no')))

    n_points <- final_table %>%
        group_by(type) %>%
        summarise(total_bigger=sum(is_1.2_order_magnitude=='bigger'),
                  total_smaller=sum(is_1.2_order_magnitude=='smaller')) %>%
        mutate(label_bigger=paste0('n = ', total_bigger), x=-Inf, y=Inf) %>%
        mutate(label_smaller=paste0('\nn = ', total_smaller), x=-Inf, y=Inf) %>%
        ungroup()
    
    return(list(final_table=final_table, n_points=n_points))
    
}

make_plot <- function(final_table, n_points) {
    
    ggplot(final_table, aes(x=PValue, y=PValue_cov)) +
        geom_point(aes(colour=is_1.2_order_magnitude), alpha=0.3, size=p.point_size_baby) +
        scale_colour_manual(values=c(p.negative_color, 'black', p.positive_color)) +
        geom_abline(color='coral3', slope=1, intercept=0) +
        geom_text(aes(x=x, y=y, label=label_bigger), data=n_points, hjust=-1, vjust=1, colour=p.negative_color, fontface='bold', size=p.label_size) +
        geom_text(aes(x=x, y=y, label=label_smaller), data=n_points, hjust=-1, vjust=1, colour=p.positive_color, fontface='bold', size=p.label_size) +
        facet_grid(.~type) +
        xlab('log10 p-value\n(original)')
        ylab('log10 p-value\n(transcriptome diversity controlled)')
        theme_noGrid() +
        theme(legend.position='top') +
        do.call(theme, theme_pars) 
    
}

differential_arid <- make_differential_table(differential$genotype)
differential_fulvestrant <- make_differential_table(differential$fulvestrant_wt)
differential_jq1 <- make_differential_table(differential$jq1_wt)
differential_tamoxifen <- make_differential_table(differential$tamoxifen_wt)

p_arid <- make_plot(differential_arid$final_table, differential_arid$n_points)
p_fulvestrant <- make_plot(differential_fulvestrant$final_table, differential_fulvestrant$n_points)
p_jq1 <- make_plot(differential_jq1$final_table, differential_jq1$n_points)
p_tamoxifen <- make_plot(differential_tamoxifen$final_table, differential_tamoxifen$n_points)

ggsave(p_file_arid, p_arid, width=5, height=2.7)
ggsave(p_file_fulvestrant, p_fulvestrant, width=5, height=2.7)
ggsave(p_file_jq1, p_jq1, width=5, height=2.7)
ggsave(p_file_tamoxifen, p_tamoxifen, width=5, height=2.7)
