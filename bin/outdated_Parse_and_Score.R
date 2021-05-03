library(dplyr)
library(tidyr)
library(readr)
library(splitstackshape)
library(stringr)

#Read-in input files
args <- commandArgs(trailingOnly = T)

#Read-in Extracted Table Files
BCSQ <- read.table(args[1], quote="\"", comment.char="")


#Parse Files
parse_BCSQ <- function(df){
  parsed <- df %>% dplyr::rename("CHROM" = "V1", "POS" = "V2", "ANNOTATION" = "V3") %>% #Separate CHROM, POS, ANNOTATION
    splitstackshape::cSplit("ANNOTATION", sep = ",", direction = "wide") #Separate Transcripts
  renamecols <- function(x){colnames(x) <- gsub("ANNOTATION", "SCRIPT", colnames(x)); x}
  pivot <- renamecols(parsed) %>%
    tidyr::pivot_longer(cols = starts_with("SCRIPT"), names_to = "TRANSCRIPT", values_to = "ANNOTATION") %>% na.omit() %>%
    tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|")%>%
    as_tibble()
  
  return(pivot)
}

parsed_BCSQ <- parse_BCSQ(BCSQ)


#Clean Table

BCSQ_con_clean <- dplyr::mutate(parsed_BCSQ, "CONSEQUENCE" = sub("*", "", parsed_BCSQ$CONSEQUENCE, fixed = TRUE))

BCSQ_p_translate <- function(df){ #Fails for AA changes that go more than one position 
  translated <- df %>%  
    dplyr::mutate("AA" = stringr::str_extract(df$AMINO_ACID_CHANGE, "[A-Z]"))%>% #Grabs the first AA 
    dplyr::mutate("ALT_AA" = stringr::str_sub(df$AMINO_ACID_CHANGE,-1)) %>% #Grabs the ALT AA
    dplyr::mutate("AA_POS" = stringr::str_extract(df$AMINO_ACID_CHANGE, "([0-9])+")) %>% #Grabs the AA posistion
    tidyr::unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|") #Unites AA, ALT
  return(translated)
} #Translate BCSQ Amino Acids 

clean_BCSQ <- BCSQ_p_translate(BCSQ_con_clean)

#Adding BLOSUM and Grantham Scores
AA_Scores <- readr::read_tsv(args[2]) #read in Amino Acid scores

clean_BCSQ <- dplyr::left_join(clean_BCSQ, AA_Scores, by = "REF_ALT_AA")


#Adding Protein Length and Calculating Percent Protein


gff_AA_Length <- readr::read_tsv(args[3]) #read in transcript AA lengths from the gff file

clean_BCSQ <- dplyr::left_join(clean_BCSQ, gff_AA_Length, by = "TRANSCRIPT")

scored_BCSQ <- clean_BCSQ %>%
  dplyr::mutate("PerProtein" = (as.numeric(AA_POS)/as.numeric(AA_Length))*100) %>%
  dplyr::mutate_if(is.numeric, round, digits=2) #Round percent protein to 2 decimal places.

#Creating BED File to Annotate VCF
prebed_BCSQ <- scored_BCSQ %>%
  dplyr::mutate("chromStart"=POS-1)

bcsq_bed <- prebed_BCSQ %>%
  dplyr::select(CHROM,chromStart,POS,BSCORE,GSCORE,PerProtein)

readr::write_tsv(bcsq_bed, "BCSQ_bed.bed", col_names = FALSE) #Header not needed
