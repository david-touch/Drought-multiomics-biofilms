

###Visualization of differentially abundant MAGs

#Load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(vegan)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data
Query <- read.delim("Data/Query.tsv", sep="\t")
MAGs_taxo <- read.csv("Data/MAGs_taxo_HQeuk.csv") %>% 
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum))
eggnog_results_bac <- read.delim("Data/all_annotations_bac.tsv", sep="\t")
eggnog_results_euk <- read.delim("Data/all_annotations_euk.tsv", sep="\t")
eggnog_results <- rbind(eggnog_results_bac, eggnog_results_euk)
rm(eggnog_results_bac, eggnog_results_euk)
DEMAGs <- read.csv("Results/DEMAGs_sign.csv") %>% dplyr::rename(Genome = MAG) %>% 
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum),
         Period = case_when(Drought == "D1" ~ 1,
                            Drought == "D2" ~ 2,
                            Drought == "D3" ~ 3,
                            TRUE ~ NA),
         
         Drought = case_when(
           log2FoldChange < 0 & Drought == "D1" ~ "B1",
           log2FoldChange > 0 & Drought == "D1" ~ "A1",
           log2FoldChange < 0 & Drought == "D2" ~ "B2",
           log2FoldChange > 0 & Drought == "D2" ~ "A2",
           log2FoldChange < 0 & Drought == "D3" ~ "B3",
           log2FoldChange > 0 & Drought == "D3" ~ "A3",
           TRUE ~ Drought))


#All droughts
RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample, levels = c("RNA1","RNA2", "RNA3", "RNA4", "RNA5", "RNA6", "RNA7","RNA8", "RNA9", "RNA10", "RNA11", "RNA12", "RNA13","RNA14", "RNA15", "RNA16", "RNA17", "RNA18"))
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought, levels = c("B1","A1", "B2", "A2", "B3", "A3"))
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

files <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata$Sample, "abundance.tsv")
names(files) <- paste0("RNA", 1:18)

tpm_list <- imap(files, ~ {
  read_tsv(.x, show_col_types = FALSE) %>%
    select(target_id, !!.y := tpm)  # Rename 'tpm' column to sample name
})
tpm_combined <- purrr::reduce(tpm_list, full_join, by = "target_id")

tpm_filtered <- tpm_combined[rowSums(tpm_combined >= 0.1) >= 2, ]

#Heatmap MAGs
tpm_df_taxo_all <- tpm_filtered %>%
  pivot_longer(cols = -target_id, names_to = "Sample", values_to = "TPM") %>%
  left_join(RNA_metadata) %>% 
  select(-Name) %>% 
  left_join(Query, by = join_by(target_id == Query)) %>% left_join(MAGs_taxo, by = join_by(Genome == MAGs)) %>% 
  group_by(Sample, Genome) %>% 
  mutate(SumTPMsample = sum(TPM)) %>% 
  ungroup() %>% 
  distinct(Sample, Genome, .keep_all = T) %>% 
  group_by(Period, Genome)%>%
  mutate(NormTPM = SumTPMsample/sum(SumTPMsample))%>%
  ungroup() %>% 
  group_by(Drought, Genome) %>% 
  mutate(AverageTPM = mean(NormTPM)) %>% 
  ungroup()

tpm_mat_taxo_all = tpm_df_taxo_all %>%
  select(Genome, NormTPM, Sample)%>%
  pivot_wider(names_from = Genome, values_from = NormTPM)%>%
  column_to_rownames("Sample")

dist_taxo_all = vegdist(t(as.matrix(tpm_mat_taxo_all)), method = "euclidian")
clus_taxo_all = hclust(dist_taxo_all, "aver")
ord_taxo_all = clus_taxo_all$order

tpm_df_taxo_all$Genome = factor(tpm_df_taxo_all$Genome, levels = colnames(tpm_mat_taxo_all)[ord_taxo_all])
tpm_df_taxo_all$Sample <- factor(tpm_df_taxo_all$Sample, levels = c("RNA1","RNA2", "RNA3", "RNA4", "RNA5", "RNA6", "RNA7","RNA8", "RNA9", "RNA10", "RNA11", "RNA12", "RNA13","RNA14", "RNA15", "RNA16", "RNA17", "RNA18"))
tpm_df_taxo_all$Phylum = factor(tpm_df_taxo_all$Phylum, levels=c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota"))
tpm_df_taxo_all$Drought = factor(tpm_df_taxo_all$Drought, levels = c("B1", "A1", "B2", "A2", "B3", "A3"))
DEMAGs$Drought = factor(DEMAGs$Drought, levels = c("B1", "A1", "B2", "A2", "B3", "A3"))
DEMAGs$Phylum = factor(DEMAGs$Phylum, levels = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota"))

ColorPhyla <- c(
  Gammaproteobacteria = "#A6CEE3",
  Alphaproteobacteria = "#1F78B4",
  Bacteroidota = "#FFFF99",
  Cyanobacteriota = "#B2DF8A", 
  Stramenopiles = "#33A02C",
  Patescibacteria = "#FB9A99",
  Actinomycetota = "#E31A1C",
  Verrucomicrobiota = "#FDBF6F", 
  Planctomycetota = "#FF7F00",  
  Chloroflexota = "#CAB2D6",
  Spirochaetota = "#6A3D9A", 
  Bdellovibrionota = "#1B9E77",
  Fibrobacterota = "#B15928", 
  Obazoa = "#1ff8ff",
  Myxococcota = "#7570B3",
  Deinococcota = "#D95f02")

Label_droughts = as_labeller(c("1" = "Drought 1 (6h)", "2" = "Drought 2 (24h)", "3" = "Drought 3 (24h)"))
Label_samples <- c("B1" = "Before", "A1" = "After",
                   "B2" = "Before", "A2" = "After",
                   "B3" = "Before", "A3" = "After")

FigureS11 <- tpm_df_taxo_all %>% 
  droplevels() %>%
  ggplot(aes(x = Drought, y = Genome)) +
  geom_tile(aes(fill= Phylum, alpha = (AverageTPM)), width = 0.95, height = 1, color = "gray", linewidth = 0.01) +
  geom_text(data = DEMAGs, aes(x = Drought, y = Genome, label = "*"), color = "red", size = 8, vjust = 0.8) +
  theme_classic() +
  scale_fill_manual(values = ColorPhyla) +
  facet_grid(Phylum ~ Period, scales = "free", space = "free", labeller = Label_droughts) +
  labs(x = NULL, y = NULL, alpha = "Normalized\ncoverage") +
  scale_x_discrete(labels = Label_samples) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"))
FigureS11

ggsave("Figures/FigureS11.pdf", FigureS11, width = 9.1, height = 12)
