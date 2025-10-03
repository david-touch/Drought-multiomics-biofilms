

###Protein extraction and identification analysis

#load libraries
library(tidyverse)
library(DEP)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data and add the Potential contaminant column
data <- read.csv('Data/combined_protein_0.3.tsv', sep='\t') 
PROT_metadata <- read.csv('Data/Protein_metadata.csv')
#annot <- read.csv("MetaG/Eggnog/HQ_eggnog.tsv", sep='\t')

eggnog_results_bac <- read.delim("Data/all_annotations_bac.tsv", sep="\t")
eggnog_results_euk <- read.delim("Data/all_annotations_euk.tsv", sep="\t")
annot <- rbind(eggnog_results_bac, eggnog_results_euk)
rm(eggnog_results_bac, eggnog_results_euk)
kegg_ko_prep <- read.csv("MetaG/kegg_ko_prep.csv")
MAGs_taxo <- read.csv("Data/MAGs_taxo_HQeuk.csv")

data <- data %>% 
  mutate(Potential.contaminant = if_else(str_detect(Protein.IDs, "CON_"), "Yes", "No")) %>% 
  filter(Reverse != "+", Potential.contaminant != "Yes")

names(data) <- names(data) %>%
  str_replace("^LFQ.Intensity\\.DT_250624_1305_Lumos_", "")

# Load annotation and add PFAMs as gene names
Prot_df <- data %>%
  separate_rows(Protein.IDs, sep = ";") %>% 
  mutate(Protein.IDs = str_trim(Protein.IDs)) %>% 
  left_join(annot, by = c("Protein.IDs" = "query")) %>% 
  rowwise() %>%
  filter(sum(c_across(starts_with("Sample"))) > 0) %>%
  ungroup() %>% 
  mutate(Common = if_else(if_all(starts_with("Sample"), ~ .x > 0), "yes", "no")) %>% 
  left_join(MAGs_taxo) %>% 
  filter(!is.na(MAGs)) %>% 
  select(Protein.IDs, Description, PFAMs, KEGG_ko, Common, MAGs, Phylum, starts_with("Sample")) 

Unique_protein <- n_distinct(Prot_df$PFAMs) #22 (without NA)
Unique_MAGs <- n_distinct(Prot_df$MAGs) #8 (without NA)
Common_protein <- n_distinct(filter(Prot_df, Common == "yes")$PFAMs) #7 

write.csv(Prot_df, "Results/Protein_results.csv", row.names = FALSE)
