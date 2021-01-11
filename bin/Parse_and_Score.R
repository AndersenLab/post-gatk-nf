#Libraries
library(tidyverse)
library(splitstackshape)
library(stringr)
library(rebus)
library(data.table)

#Read-in Extracted Table Files
BCSQ <- read.table("BCSQ.vcf", quote="\"", comment.char="")


#Parse Files
parse_BCSQ <- function(df){
  parsed <- df %>% rename("CHROM" = "V1", "POS" = "V2", "ANNOTATION" = "V3") %>% #Separate CHROM, POS, ANNOTATION
    cSplit("ANNOTATION", sep = ",", direction = "wide") #Separate Transcripts
  renamecols <- function(x){colnames(x) <- gsub("ANNOTATION", "SCRIPT", colnames(x)); x}
  pivot <- renamecols(parsed) %>%
    pivot_longer(cols = starts_with("SCRIPT"), names_to = "TRANSCRIPT", values_to = "ANNOTATION") %>% na.omit() %>%
    separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|")%>%
    as_tibble()
  
  return(pivot)
}

parsed_BCSQ <- parse_BCSQ(BCSQ)


#Clean Table

BCSQ_con_clean <- mutate(parsed_BCSQ, "CONSEQUENCE" = sub("*", "", parsed_BCSQ$CONSEQUENCE, fixed = TRUE))

BCSQ_p_translate <- function(df){ #Fails for AA changes that go more than one position 
  translated <- df %>%  
    mutate("AA" = str_extract(df$AMINO_ACID_CHANGE, "[A-Z]"))%>% #Grabs the first AA 
    mutate("ALT_AA" = str_sub(df$AMINO_ACID_CHANGE,-1)) %>% #Grabs the ALT AA
    mutate("AA_POS" = str_extract(df$AMINO_ACID_CHANGE, "([0-9])+")) %>% #Grabs the AA posistion
    unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|") #Unites AA, ALT
  return(translated)
} #Translate BCSQ Amino Acids 

clean_BCSQ <- BCSQ_p_translate(BCSQ_con_clean)


#Adding BLOSUM and Grantham Scores

AA_Scores <- AA_Scores <- read_tsv("AA_Scores.tsv") #read in Amino Acid scores

clean_BCSQ <- left_join(clean_BCSQ, AA_Scores, by = "REF_ALT_AA")


#Adding Protein Length and Calculating Percent Protein
gff_AA_Length <- read_tsv("gff_AA_Length.tsv") #read in transcript AA lengths from the gff file

clean_BCSQ <- left_join(clean_BCSQ, gff_AA_Length, by = "TRANSCRIPT")

scored_BCSQ <- clean_BCSQ %>%
  mutate("PerProtein" = (AA_POS/AA_Length)*100) %>%
  mutate_if(is.numeric, round, digits=2) #Round percent protein to 2 decimal places.

#Creating BED File to Annotate VCF
prebed_BCSQ <- scored_BCSQ %>%
  mutate("chromStart"=POS-1)

bcsq_bed <- prebed_BCSQ %>%
  select(CHROM,chromStart,POS,BSCORE,GSCORE,PerProtein)

write_tsv(bcsq_bed, "BCSQ_bed.bed", col_names = FALSE) #Header not needed
