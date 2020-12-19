GTEX_TISSUES_LOW_SAMPLES <- c('Kidney_Medulla', 'Fallopian_Tube', 'Cervix_Endocervix', 'Cervix_Ectocervix', 'Bladder')

correct_tissue_names <- function(tissues) {
    
    tissues <- gsub('gtex_', '', tissues)
    
    convert_tis <- setNames(nm = c("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Small_Intestine_Terminal_Ileum", 
                                   "Muscle_Skeletal", "Breast_Mammary_Tissue", "Whole_Blood", "Cells_EBVtransformed_lymphocytes"),
                            object = c("Skin sun exposed", "Skin sun protected", "Small Intestine", 
                                       "Muscle", "Breast", "Blood", "Cells EBV-transformed lymphocytes"))
    
    # renaming
    for(i in seq_along(convert_tis)) {
        tissues[tissues == names(convert_tis[i])] <- convert_tis[tissues[tissues == names(convert_tis[i])]]
    }
    
    # Eliminating under score
    tissues <- gsub("_", " ", tissues)
    
    return(tissues)
    
}
