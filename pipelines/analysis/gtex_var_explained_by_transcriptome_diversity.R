source('../../R/misc.R')
library('dplyr')
path <- '/scratch/users/paedugar/transcriptome_diversity/var_explained_by_transcriptome_diversity/pca_method/from_raw/'

all_files <- list.files(path, full.names='T', recursive=T, pattern='var_explained')
all_files <- all_files[grep('gtex', all_files)]

x <- concatenate_table_files2(all_files, id_names=all_files)
x$method <- ifelse(grepl('tpm', x$id_names), 'tpm', 'tmm')
x$id_names <- basename(dirname(x$id_names))

# Median per method
med <- x %>%
    group_by(method) %>%
    summarise(med=median(var_explained)) %>%
    ungroup()

print(med)

x <- x %>% 
    group_by(method) %>%
    arrange(-var_explained) %>%
    ungroup()

print(as.data.frame(x))
