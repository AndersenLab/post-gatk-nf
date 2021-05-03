library(tidyverse)

join_masks <- function(mask_file = NULL){
    cluster <- NA
    cluster_start <- NA
    
    for (i in c(1:nrow(mask_file)-1)) {
        
        if(i==1){
            
            cluster[i] <- ifelse(mask_file$END_BIN[i] == mask_file$START_BIN[i+1], 'yes', 'no')
            cluster_start[i] <- mask_file$START_BIN[i]
            
        } else {
            
            cluster[i] <- ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i] | mask_file$END_BIN[i] == mask_file$START_BIN[i+1], 'yes', 'no')
            cluster_start[i] <- ifelse(mask_file$CHROM[i] == mask_file$CHROM[i-1] & cluster[i] == 'yes',
                                       ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i] | mask_file$END_BIN[i] == mask_file$START_BIN[i+1],
                                              ifelse(cluster[i-1] == "no", mask_file$START_BIN[i],
                                                     ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i], cluster_start[i-1], mask_file$START_BIN[i])), mask_file$START_BIN[i]), mask_file$START_BIN[i])
            
        }
    }
    
    cluster[nrow(mask_file)] <- NA
    cluster_start[nrow(mask_file)] <- NA
    
    mask_file_cluster <- data.frame(mask_file, cluster, cluster_start) %>%
        dplyr::group_by(CHROM, cluster_start) %>%
        dplyr::mutate(cluster_size=ifelse(is.na(cluster_start), 1000, n()*1000)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(cluster_end=cluster_start+cluster_size)
    
    
}

masked_directory <- "."
#mask_files <- list.files(masked_directory)

filename_list = list.files(path = masked_directory, pattern="*_Mask_DF.tsv", full.names = T)
# path do not need the final "/". full.names add path to the list

file_list <- lapply(filename_list, read_delim, "\t", escape_double = FALSE, trim_ws = TRUE)

names(file_list) <- stringr::str_replace(basename(filename_list), pattern = "_Mask_DF.tsv", replacement = "")

population_all = lapply(names(file_list), function(nm) mutate(file_list[[nm]], STRAIN=nm)) %>% purrr::reduce(bind_rows)

strain_count <- length(filename_list)

write.table(population_all, file = "divergent_output_all_strains_all_bins.tsv")