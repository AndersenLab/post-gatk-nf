print("Starting second script")

library(tidyverse)
masked_directory <- "."
filename_list = list.files(path = masked_directory, pattern="*_Mask_DF.tsv", full.names = T)
strain_count <- length(filename_list)

population_all <- read.table(file = "divergent_output_all_strains_all_bins.tsv")

print("finished reading in data")

# COVERAGE
df_coverage <- population_all %>%
    dplyr::group_by(STRAIN) %>%
    dplyr::summarise(mean_cov_genome = mean(COVERAGE)) %>%
    dplyr::ungroup()

df_coverage_fraction <- population_all %>%
    dplyr::left_join(., df_coverage, by = "STRAIN") %>%
    dplyr::mutate(fraction_cov = COVERAGE/mean_cov_genome)

print("Finished step 1")

## long-read based optimized paramters (2020/05/11)

cov=0.35
ct=15
cluster_threshold=9000
gap_cluster_threshold=9000
gap=5000

population_all_div <- df_coverage_fraction  %>%
    dplyr::mutate(window_ID = paste(CHROM, START_BIN, END_BIN)) %>%
    dplyr::mutate(div_count = ifelse(VAR_COUNT>ct, "Count", "Pass"), div_lowcov = ifelse(fraction_cov<cov, "LowCoverage", "Pass")) %>%
    dplyr::mutate(div_class = ifelse(paste(div_count, div_lowcov, sep="_") == "Pass_Pass", "Pass", "Divergent")) %>%
    dplyr::mutate(div_twoflank = ifelse((grepl("Pass", dplyr::lag(div_class)) | grepl("Pass", dplyr::lead(div_class))) | !grepl("Pass", div_class), "Pass", "TwoFlank")) %>%
    dplyr::mutate(div_id = paste(ifelse(div_count=="Pass","",div_count), ifelse(div_lowcov=="Pass","",div_lowcov), ifelse(div_twoflank=="Pass","",div_twoflank), sep=""), div_class = ifelse(paste(div_count, div_lowcov, div_twoflank, sep="_") == "Pass_Pass_Pass", "Pass", "Divergent"))

print("Finished step 2")

# JOINING
cluster_start <- NA
cluster_end <- NA

df_div_cluster <- population_all_div %>%
    dplyr::group_by(STRAIN, CHROM) %>%
    dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
    dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
    dplyr::ungroup() %>%
    dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)

df_div_cluster_filtered <- df_div_cluster %>%
    dplyr::filter(cluster == T)

print("finished step 3")

cluster_start_vec <- df_div_cluster_filtered$cluster_start
start_bin_vec <- df_div_cluster_filtered$START_BIN

for (i in 1:length(cluster_start_vec)) {
    cluster_start_vec[i] <- ifelse(start_bin_vec[i]==0, 0, 
                                   ifelse(!is.na(cluster_start_vec[i]), cluster_start_vec[i], cluster_start_vec[i-1]))
    
}

print("finished step 4")

df_div_cluster_filtered <- df_div_cluster_filtered %>%
    dplyr::select(-cluster_start) %>%
    dplyr::mutate(cluster_start=cluster_start_vec)

df_div_cluster_size <- df_div_cluster %>%
    dplyr::select(-cluster_start) %>%
    dplyr::left_join(., dplyr::select(df_div_cluster_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
    dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
    dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_end=cluster_start+cluster_size)

print("finished step 5")

# FILTERING AND GAP FILLING
divergent_gap <- function(df=df_joined_masks_freq, threshold = gap_cluster_threshold, gap_threshold=gap) {
    
    df_gap <- df %>%
        dplyr::filter(cluster_size >= threshold) %>%
        dplyr::distinct(STRAIN, CHROM, cluster_start,cluster_size, cluster_end)
    
    gap <- NA
    gap_size <- NA
    gap_join <- NA
    extended_cluster_end <- NA
    extended_cluster_size <-NA
    
    for (i in c(1:nrow(df_gap)-1)) {
        
        gap[i] <- ifelse(df_gap$STRAIN[i] == df_gap$STRAIN[i+1] & df_gap$CHROM[i] == df_gap$CHROM[i+1], 'yes', NA)
        gap_size[i] <- ifelse(gap[i]=="yes", df_gap$cluster_start[i+1]-df_gap$cluster_end[i], NA)
        gap_join[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, 'join', NA)
        extended_cluster_end[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, df_gap$cluster_end[i+1],df_gap$cluster_end[i])
        extended_cluster_size[i] <- ifelse(gap[i]=="yes" & gap_size[i] <= gap_threshold & gap_size[i] > 0, extended_cluster_end[i]-df_gap$cluster_start[i], df_gap$cluster_size[i])
        
    }
    
    gap[nrow(df_gap)] <- NA
    gap_size[nrow(df_gap)] <- NA
    gap_join[nrow(df_gap)] <- NA
    extended_cluster_end[nrow(df_gap)] <- NA
    extended_cluster_size[nrow(df_gap)] <- NA
    
    df_gap_size <- data.frame(df_gap, gap, gap_size, gap_join, extended_cluster_end, extended_cluster_size) %>%
        dplyr::mutate(threshold = threshold, threshold_pass = ifelse(cluster_size >= threshold, T, F)) %>%
        dplyr::group_by(STRAIN, CHROM, extended_cluster_end) %>%
        dplyr::mutate(extended_cluster_start = min(cluster_start)) %>%
        dplyr::mutate(extended_cluster_size = extended_cluster_end - extended_cluster_start) %>%
        dplyr::ungroup()
    
    df_gap_size_join <- df %>%
        dplyr::left_join(., df_gap_size, by=c('STRAIN', 'CHROM', 'cluster_start', 'cluster_size', 'cluster_end'))
    
    
    return(list(df_gap_size, df_gap_size_join))
    
}

df_gap_bind_list <- divergent_gap(df=df_div_cluster_size, threshold = gap_cluster_threshold, gap_threshold=gap)

df_gap_bind <- df_gap_bind_list[[1]]
df_gap_bind_all_bins <- df_gap_bind_list[[2]] %>%
    dplyr::mutate(div_id = ifelse(div_class == "Divergent", div_id, 
                                  ifelse(dplyr::lag(END_BIN, 1) == dplyr::lag(cluster_end, 1) & dplyr::lag(gap_join, 1) == 'join' | 
                                             dplyr::lag(END_BIN, 2) == dplyr::lag(cluster_end, 2) & dplyr::lag(gap_join, 2) == "join" | 
                                             dplyr::lag(END_BIN, 3) == dplyr::lag(cluster_end, 3) & dplyr::lag(gap_join, 3) == "join" |
                                             dplyr::lag(END_BIN, 4) == dplyr::lag(cluster_end, 4) & dplyr::lag(gap_join, 4) == "join", "gap", div_id))) %>%
    dplyr::mutate(div_class = ifelse(div_id == "gap", "Divergent", div_class)) %>%
    dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id)

print("finished step 6")

# RECLUSTERING
df_div_cluster_gapjoined <- df_gap_bind_all_bins %>%
    dplyr::group_by(STRAIN, CHROM) %>%
    dplyr::mutate(cluster = ifelse(((grepl("Divergent", dplyr::lag(div_class)) | grepl("Divergent", dplyr::lead(div_class)))) & grepl("Divergent", div_class), T, F)) %>%
    dplyr::mutate(cluster_start = ifelse(dplyr::lag(cluster) == F & cluster == T & dplyr::lead(cluster) == T, START_BIN,NA)) %>%
    dplyr::ungroup() %>%
    dplyr::select(STRAIN, CHROM, START_BIN, END_BIN, window_ID, div_class, div_id, cluster, cluster_start)

df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined %>%
    dplyr::filter(cluster == T)

print("finished step 7")

cluster_start_gapjoined_vec <- df_div_cluster_gapjoined_filtered$cluster_start
start_bin_gapjoined_vec <- df_div_cluster_gapjoined_filtered$START_BIN

for (i in 1:length(cluster_start_gapjoined_vec)) {
    cluster_start_gapjoined_vec[i] <- ifelse(start_bin_gapjoined_vec[i]==0, 0, 
                                             ifelse(!is.na(cluster_start_gapjoined_vec[i]), cluster_start_gapjoined_vec[i], cluster_start_gapjoined_vec[i-1]))
}

df_div_cluster_gapjoined_filtered <- df_div_cluster_gapjoined_filtered %>%
    dplyr::select(-cluster_start) %>%
    dplyr::mutate(cluster_start=cluster_start_gapjoined_vec)

print("finished step 8")

df_div_cluster_gapjoined_size <- df_div_cluster_gapjoined %>%
    dplyr::select(-cluster_start) %>%
    dplyr::left_join(., dplyr::select(df_div_cluster_gapjoined_filtered, STRAIN, window_ID, cluster_start), by=c('STRAIN','window_ID')) %>%
    dplyr::group_by(STRAIN, CHROM, cluster_start) %>%
    dplyr::mutate(cluster_size=ifelse(div_class %in% c("Pass", "N2_fp"), NA, ifelse(is.na(cluster_start), 1000, n()*1000))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_end=cluster_start+cluster_size)

df_divergent_final <- df_div_cluster_gapjoined_size %>%
    dplyr::filter(cluster_size >= cluster_threshold) %>%
    dplyr::group_by(window_ID) %>%
    dplyr::mutate(freq=n()/strain_count, mwf = ifelse(freq > 0.5, 1-freq, freq)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(freq_bin=ifelse(freq < 0.01, "Rare", ifelse(freq < 0.05, "Intermediate", "Common"))) %>%
    dplyr::mutate(is_count = ifelse(div_id %in% c("Count", "CountLowCoverage"), 1, 0)) %>%
    dplyr::group_by(STRAIN, CHROM, cluster_start, cluster_end) %>%
    dplyr::mutate(flag_del = ifelse(max(is_count)==0, "Del", "Pass")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(flag_del == "Pass")

print("finished divergent final")

write.table(df_divergent_final, file = "df_divergent_final.tsv")

# SUMMARY PER ISOTYPE
df_divergent_final_isotype <- df_divergent_final %>%
    dplyr::distinct(STRAIN, CHROM, cluster_start, cluster_size, cluster_end)

print("finished isotypes")

write.csv(df_divergent_final_isotype, file="DataS3_divergent_regions_isotypes.csv")

df_divergent_final_isotype %>% select(CHROM, cluster_start, cluster_end, STRAIN) %>% arrange(CHROM, cluster_start) %>% write.table("divergent_regions_strain.bed", quote=F, col.names = F, row.names=F, sep="\t")