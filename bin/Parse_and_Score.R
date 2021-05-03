#Libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#Read-in input files
args <- commandArgs(trailingOnly = T)

#Read-in Extracted Table Files
BCSQ <- read.table(args[1], quote="\"", comment.char="", stringsAsFactors = FALSE)


#Parse Files
parse_VCF <- function(df){
  parsed <- df %>% 
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "ANNOTATION" = "V3") %>%
    tidyr::separate_rows(ANNOTATION, sep=",")%>%
    tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") 
}

parsed_BCSQ <- parse_VCF(BCSQ)


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

write.csv(clean_BCSQ, "clean_BCSQ.csv") #Write out clean BCSQ

#Adding BLOSUM and Grantham Scores

clean_BCSQ <- readr::read_csv("clean_BCSQ.csv") #Read in clean BCSQ

AA_Scores <- readr::read_tsv(args[2]) #read in Amino Acid scores

clean_BCSQ <- dplyr::left_join(clean_BCSQ, AA_Scores, by = "REF_ALT_AA")


#Adding Protein Length and Calculating Percent Protein


gff_AA_Length <- readr::read_tsv(args[3]) #read in transcript AA lengths from the gff file

clean_BCSQ <- dplyr::left_join(clean_BCSQ, gff_AA_Length, by = "TRANSCRIPT")

scored_BCSQ <- clean_BCSQ %>%
  dplyr::mutate("PerProtein" = (AA_POS/AA_Length)*100) %>%
  dplyr::mutate_if(is.numeric, round, digits=2) #Round percent protein to 2 decimal places.

write.csv(scored_BCSQ, "scored_BCSQ.csv") #Write out for analysis

#Creating BED File to Annotate VCF

scored_BCSQ <- readr::read_csv("scored_BCSQ.csv")

prebed_BCSQ <- scored_BCSQ %>%
  dplyr::mutate("chromStart"=POS-1)

bcsq_bed <- prebed_BCSQ %>%
  dplyr::select(CHROM,chromStart,POS,BSCORE,GSCORE,PerProtein)

readr::write_tsv(bcsq_bed, "BCSQ_bed.bed", col_names = FALSE) #Header not needed
