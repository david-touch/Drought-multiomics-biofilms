

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
Drought_gene_sign <- read.csv("Results/Drought_df_sign.csv")
kegg_ko_prep <- read.csv("Data/kegg_ko_prep.csv")


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

#Counts per category
DEG_func_filt <- Drought_gene_sign %>% 
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:")) %>% 
  filter(KEGG_ko != "NA" & KEGG_ko != "-") %>% 
  left_join(kegg_ko_prep, by = join_by(KEGG_ko == Ko), relationship = "many-to-many") %>% 
  filter(Category != "Organismal Systems", 
         Category != "Human Diseases",
         Subcategory != "Poorly characterized") %>% 
  filter(!(Kingdom == "Bacteria" & Subcategory %in% c("Transport and catabolism"))) %>% 
  filter(!(Kingdom == "Bacteria" & Pathway %in% c("Ubiquitin mediated proteolysis"))) %>% 
  filter(!(Kingdom == "Bacteria" & Subcategory == "Cell growth and death" & Pathway != "Cell cycle - Caulobacter")) %>% 
  mutate(Category = case_when(Category %in% c("Brite Hierarchies", "Not Included in Pathway or Brite") ~ "Others Functions",
                              TRUE ~ Category))

DEG_func_count <- DEG_func_filt %>%
  summarise(KEGG_ko = n_distinct(KEGG_ko),
            Category = n_distinct(Category),
            Subcategory = n_distinct(Subcategory),
            Pathway = n_distinct(Pathway))

DEG_fun_df <- DEG_func_filt %>% 
  add_count(DA_Period, Category, name = "sumSubCat") %>% 
  distinct(DA_Period, Category, .keep_all = T)

DEG_cat_count <- DEG_func_filt %>% 
  add_count(Category, name = "sumCat") %>% 
  distinct(Category, .keep_all = T)

DEG_path_count <- DEG_func_filt %>% 
  add_count(Pathway, name = "sumPath") %>% 
  distinct(Pathway, .keep_all = T)

Label_droughts = as_labeller(c("D1" = "Drought 1 (6h)", "D2" = "Drought 2 (24h)", "D3" = "Drought 3 (24h)"))

barplot_DEG_func <- DEG_fun_df %>% 
  group_by(Drought) %>% 
  complete(Category, DA_Period, fill = list(sumSubCat = 1)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Category, y = sumSubCat, fill = DA_Period)) +
  geom_col(position = "dodge") +
  labs(
    x = "KEGG category",
    y = "Number of DEGs",
    fill = "DEGs",
  ) +
  scale_fill_manual(values = c(
    "Down-regulated D1" = "#a6cee3",
    "Down-regulated D2" = "#1f78b4",
    "Down-regulated D3" = "#295184",
    "Up-regulated D1" = "#f1948a",
    "Up-regulated D2" = "#ec7063",
    "Up-regulated D3" = "#b03a2e"),
    breaks = c("Down-regulated D2", "Up-regulated D2"),
    labels = c("Pre-drought", "Post-drought")) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(~ Drought, labeller = Label_droughts) + 
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y =  element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black", size = 0.5, linetype = "solid"),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 16, face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray80", size = 0.2))
barplot_DEG_func

ggsave("Figures/FigureS12.pdf", barplot_DEG_func, width = 18.2, height = 9.1)

#Heatmap function
tpm_df_func <- tpm_filtered %>%
  filter(target_id %in% Drought_gene_sign$Gene) %>% 
  pivot_longer(cols = -target_id, names_to = "Sample", values_to = "TPM") %>%
  left_join(RNA_metadata) %>% 
  select(-Name) %>% 
  left_join(eggnog_results, by = join_by(target_id == query)) %>%
  select(target_id, Sample, TPM, Time, Drought, State, Period, Succession, MAGs, KEGG_ko) %>% 
  left_join(MAGs_taxo) %>% 
  filter(KEGG_ko != "-") %>% 
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:")) %>% 
  left_join(kegg_ko_prep, by = join_by(KEGG_ko == Ko), relationship = "many-to-many")


tpm_df_func_filt <- tpm_df_func %>% 
  filter(Category != "Organismal Systems", 
         Category != "Human Diseases",
         Subcategory != "Poorly characterized") %>% 
  filter(!(Kingdom == "Bacteria" & Subcategory %in% c("Transport and catabolism"))) %>% 
  filter(!(Kingdom == "Bacteria" & Pathway %in% c("Ubiquitin mediated proteolysis"))) %>% 
  filter(!(Kingdom == "Bacteria" & Subcategory == "Cell growth and death" & Pathway != "Cell cycle - Caulobacter")) %>% 
  mutate(Category = case_when(Category %in% c("Brite Hierarchies", "Not Included in Pathway or Brite") ~ "Others Functions",
                              TRUE ~ Category))

#By subcategory
tpm_df_func_filt_norm_sub <- tpm_df_func_filt %>% 
  group_by(Sample, Subcategory) %>% 
  mutate(SumTPMsample = sum(TPM)) %>%
  ungroup() %>% 
  distinct(Sample, Subcategory, .keep_all = T) %>% 
  group_by(Period, Subcategory) %>% 
  mutate(NormTPM = SumTPMsample/sum(SumTPMsample)) %>%
  ungroup() %>% 
  group_by(Drought, Subcategory) %>% 
  mutate(AverageTPM = mean(NormTPM)) %>% 
  ungroup() %>% 
  mutate(AverageTPM = ifelse(is.nan(AverageTPM), 0, AverageTPM),
         NormTPM = ifelse(is.nan(NormTPM), 0, NormTPM))

tpm_mat_func_sub = tpm_df_func_filt_norm_sub %>%
  select(Subcategory, NormTPM, Sample)%>%
  pivot_wider(names_from = Subcategory, values_from = NormTPM)%>%
  column_to_rownames("Sample")

dist_func_sub = vegdist(t(as.matrix(tpm_mat_func_sub)), method = "euclidian")
clus_func_sub = hclust(dist_func_sub, "aver")
ord_func_sub = clus_func_sub$order

tpm_df_func_filt_norm_sub$Subcategory = factor(tpm_df_func_filt_norm_sub$Subcategory, levels = colnames(tpm_mat_func_sub)[ord_func_sub])

#By pathway
tpm_df_func_filt_norm_pat <- tpm_df_func_filt %>% 
  group_by(Sample, Pathway) %>% 
  mutate(SumTPMsample = sum(TPM)) %>%
  ungroup() %>% 
  distinct(Sample, Pathway, .keep_all = T) %>% 
  group_by(Period, Pathway) %>% 
  mutate(NormTPM = SumTPMsample/sum(SumTPMsample)) %>%
  ungroup() %>% 
  group_by(Drought, Pathway) %>% 
  mutate(AverageTPM = mean(NormTPM)) %>% 
  ungroup() %>% 
  mutate(AverageTPM = ifelse(is.nan(AverageTPM), 0, AverageTPM),
         NormTPM = ifelse(is.nan(NormTPM), 0, NormTPM))

tpm_mat_func_pat = tpm_df_func_filt_norm_pat %>%
  select(Pathway, NormTPM, Sample)%>%
  pivot_wider(names_from = Pathway, values_from = NormTPM)%>%
  column_to_rownames("Sample")

dist_func_pat = vegdist(t(as.matrix(tpm_mat_func_pat)), method = "euclidian")
clus_func_pat = hclust(dist_func_pat, "aver")
ord_func_pat = clus_func_pat$order

tpm_df_func_filt_norm_pat$Pathway = factor(tpm_df_func_filt_norm_pat$Pathway, levels = colnames(tpm_mat_func_pat)[ord_func_pat])

#Plotting
Label_droughts = as_labeller(c("1" = "Drought 1 (6h)", "2" = "Drought 2 (24h)", "3" = "Drought 3 (24h)"))
Label_samples <- c("B1" = "Before", "A1" = "After",
                   "B2" = "Before", "A2" = "After",
                   "B3" = "Before", "A3" = "After")


heatmap_func_subcat_sign <- tpm_df_func_filt_norm_sub %>% 
  droplevels() %>%
  ggplot(aes(x = State, y = Subcategory)) +
  theme_classic() +
  geom_tile(aes(fill= Category, alpha = (NormTPM)), width = 0.95, height = 1, color = "gray", linewidth = 0.01) +
  facet_grid(Category ~ Period, scales = "free", space = "free", labeller = Label_droughts) +
  labs(x = NULL, y = NULL, alpha = "Normalized\ncoverage") +
  scale_x_discrete(labels = Label_samples) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 14),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"))
heatmap_func_subcat_sign

ggsave("Figures/FigureS13.pdf", heatmap_func_subcat_sign, width = 14, height = 12)

heatmap_func_pat_sign <- tpm_df_func_filt_norm_pat %>% 
  droplevels() %>%
  ggplot(aes(x = State, y = Pathway)) +
  theme_classic() +
  geom_tile(aes(fill= Category, alpha = (NormTPM)), width = 0.95, height = 1, color = "gray", linewidth = 0.01) +
  facet_grid(Category ~ Period, scales = "free", space = "free", labeller = Label_droughts) +
  labs(x = NULL, y = NULL, alpha = "Normalized\ncoverage") +
  scale_x_discrete(labels = Label_samples) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"))
heatmap_func_pat_sign

ggsave("Figures/FigureS14.pdf", heatmap_func_pat_sign, width = 9.1, height = 12)

##Functional redundancy
FR_df <- eggnog_results %>% 
  left_join(MAGs_taxo) %>% 
  filter(Phylum %in% c("Stramenopiles", "Cyanobacteriota")) %>% 
  select(MAGs, KEGG_ko, Phylum) %>% 
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:")) %>% 
  filter(KEGG_ko != "NA" & KEGG_ko != "-") %>% 
  group_by(Phylum) %>% 
  distinct(KEGG_ko, .keep_all = T)
#File to upload in KEGG Mapper - Reconstruct: https://www.genome.jp/kegg/mapper/reconstruct.html
