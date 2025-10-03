

###Differential expression analysis on MAGs and "Global-scaled"

#Load libraries
library(tximport)
library(readr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggalluvial)
library(ggh4x)
library(ggrepel)
  

#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)

#Load data
eggnog_results_bac <- read.delim("Data/all_annotations_bac.tsv", sep="\t")
eggnog_results_euk <- read.delim("Data/all_annotations_euk.tsv", sep="\t")
eggnog_results <- rbind(eggnog_results_bac, eggnog_results_euk)
rm(eggnog_results_bac, eggnog_results_euk)
MAGs_taxo <- read.csv("Data/MAGs_taxo_HQeuk.csv")

#Metadata preparation
RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample)
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought)
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

#List negative control
file <- file.path(".", "Data", "Kallisto", "HQ_Genes", "RNANEG", "abundance.tsv")
names(file) <- paste0("RNANEG")

NEG_list <- imap(file, ~ {
  read_tsv(.x, show_col_types = FALSE) %>%
    select(target_id, !!.y := tpm)  # Rename 'tpm' column to sample name
})
NEG_list <- purrr::reduce(NEG_list, full_join, by = "target_id") %>%  filter(RNANEG > 0)


##DESeq2 analysis
RNA_metadata <- RNA_metadata %>% filter(Sample != "RNANEG")

#Filter out reads: more than 10 reads and in at least 2 samples.
files <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata$Sample, "abundance.h5")
names(files) <-c(paste0("RNA", 1:18))
txi.kallisto_all <- tximport(files, type = "kallisto", txOut = TRUE)
RNA_metadata_mod <- RNA_metadata %>% mutate(ID = Sample) %>% column_to_rownames("Sample") %>% select(Time, Drought, State, Period, ID)

dds_all <- DESeqDataSetFromTximport(txi.kallisto_all, RNA_metadata_mod, ~ State)
rm(txi.kallisto_all)
keep_all <-rowSums(counts(dds_all) >= 1) >= 2 & rowSums(counts(dds_all)) >= 10
dds_all_filt <- dds_all[keep_all,]

#Transcriptomic PCA
dds_all_res <- DESeq(dds_all_filt)
res_all <- results(dds_all_res, alpha = 0.05, , pAdjustMethod = "BH")
res_all_ordered <- res_all[order(res_all$pvalue),]
summary(res_all_ordered)
Dall <- as.data.frame(res_all_ordered) %>% rownames_to_column("Gene")
#6 up, 0 down genes, with alpha = 0.05 and filtering 

vsd_all <- vst(dds_all_res, blind=TRUE)
pcaData_all <- plotPCA(vsd_all, intgroup=c("Drought", "Period", "ID"), returnData=T)
percentVar_all <- round(100 * attr(pcaData_all, "percentVar"))

pcaData_all$Drought <- factor(pcaData_all$Drought, levels = c("B1", "A1", "B2", "A2", "B3", "A3"))

FigureS9 <- ggplot(pcaData_all, aes(PC1, PC2, color=Drought, shape=Period)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar_all[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_all[2],"% variance")) + 
  theme_classic()+
  scale_color_manual(values = c("#a6cee3","#f1948a", "#1f78b4", "#ec7063", "#295184","#b03a2e"), labels = c("Pre-drought 1 (6h)", "Post-drought 1 (6h)", "Pre-drought 2 (24h)", "Post-drought 2 (24h)", "Pre-drought 3 (24h)", "Post-drought 3 (24h)")) +
  geom_text(aes(label=ID), color = "black", vjust = -1.2) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold")) +
    guides(shape = "none",
           color = guide_legend(override.aes = list(shape = c(16,16,17,17,15,15))))
FigureS9
ggsave("Figures/FigureS9.pdf", FigureS9, width = 9.1, height = 9.1)


##DEGenes

#DEG First drought
RNA_metadata_D1 <- RNA_metadata %>% filter(Time == 4 | Time == 5)
files_D1 <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata_D1$Sample, "abundance.h5")
names(files_D1) <- paste0("RNA", 1:6)
txi.kallisto_D1 <- tximport(files_D1, type = "kallisto", txOut = TRUE)
RNA_metadata_D1_mod <- RNA_metadata_D1 %>% column_to_rownames("Sample") %>% select(Time, Drought, State)
dds_D1 <- DESeqDataSetFromTximport(txi.kallisto_D1, RNA_metadata_D1_mod, ~ State)
rm(txi.kallisto_D1)
dds_D1 <- dds_D1[!(rownames(dds_D1) %in% NEG_list$target_id), ]
keep_D1 <- rowSums(counts(dds_D1) >= 1) >= 2 & rowSums(counts(dds_D1)) >= 10
dds_D1 <- dds_D1[keep_D1, ]

dds_D1_res <- DESeq(dds_D1)
res_D1 <- results(dds_D1_res, alpha = 0.05, pAdjustMethod = "BH", contrast=c("State","After","Before"))
res_D1_ordered <- res_D1[order(res_D1$pvalue),]
summary(res_D1_ordered)
D1 <- as.data.frame(res_D1_ordered) %>% rownames_to_column("Gene") %>% mutate(Drought = "D1")
#199 up, 408 down genes, with alpha = 0.05 and filtering


#DEG Second drought
RNA_metadata_D2 <- RNA_metadata %>% filter(Time == 8 | Time == 9)
files_D2 <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata_D2$Sample, "abundance.h5")
names(files_D2) <- paste0("RNA", 7:12)
txi.kallisto_D2 <- tximport(files_D2, type = "kallisto", txOut = TRUE)
RNA_metadata_D2_mod <- RNA_metadata_D2 %>% column_to_rownames("Sample") %>% select(Time, Drought, State)
dds_D2 <- DESeqDataSetFromTximport(txi.kallisto_D2, RNA_metadata_D2, ~ State)
rm(txi.kallisto_D2)
dds_D2 <- dds_D2[!(rownames(dds_D2) %in% NEG_list$target_id), ]
keep_D2 <- rowSums(counts(dds_D2) >= 1) >= 2 & rowSums(counts(dds_D2)) >= 10
dds_D2 <- dds_D2[keep_D2, ]

dds_D2_res <- DESeq(dds_D2)
res_D2 <- results(dds_D2_res, alpha = 0.05, pAdjustMethod = "BH", contrast=c("State","After","Before"))
res_D2_ordered <- res_D2[order(res_D2$pvalue),]
summary(res_D2_ordered)
D2 <- as.data.frame(res_D2_ordered) %>% rownames_to_column("Gene") %>% mutate(Drought = "D2")
#6200 up, 2298 down genes, with alpha = 0.05 and filtering 

#DEG Third drought
RNA_metadata_D3 <- RNA_metadata %>% filter(Time == 17 | Time == 18)
files_D3 <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata_D3$Sample, "abundance.h5")
names(files_D3) <- paste0("RNA", 13:18)
txi.kallisto_D3 <- tximport(files_D3, type = "kallisto", txOut = TRUE)
RNA_metadata_D3_mod <- RNA_metadata_D3 %>% column_to_rownames("Sample") %>% select(Time, Drought, State)
dds_D3 <- DESeqDataSetFromTximport(txi.kallisto_D3, RNA_metadata_D3, ~ State)
rm(txi.kallisto_D3)
dds_D3 <- dds_D3[!(rownames(dds_D3) %in% NEG_list$target_id), ]
keep_D3 <- rowSums(counts(dds_D3) >= 1) >= 2 & rowSums(counts(dds_D3)) >= 10
dds_D3 <- dds_D3[keep_D3, ]

dds_D3_res <- DESeq(dds_D3)
res_D3 <- results(dds_D3_res, alpha = 0.05, pAdjustMethod = "BH")
res_D3_ordered <- res_D3[order(res_D3$pvalue),]
summary(res_D3_ordered)
D3 <- as.data.frame(res_D3_ordered) %>% rownames_to_column("Gene") %>% mutate(Drought = "D3")
#6 up, 12 down genes, with alpha = 0.05 and filtering 


##DEMAGs
design <- as.formula(~ State)

#DEM First drought
counts_d1 <- counts(dds_D1) %>% data.frame() %>% rownames_to_column("Query")

df_D1 <- counts_d1 %>% left_join(eggnog_results, by = join_by(Query == query)) %>% filter(!is.na(MAGs)) %>%  
  group_by(MAGs) %>% 
  summarise(across(starts_with("RNA"), sum)) %>% 
  ungroup() %>% 
  column_to_rownames("MAGs") %>% 
  round() %>% 
  as.matrix()

D1.raw <- DESeqDataSetFromMatrix(countData = df_D1,
                                 colData = RNA_metadata_D1,
                                 design = design)
D1_obj <- DESeq(D1.raw)
res_D1_t <- results(D1_obj, alpha = 0.05, pAdjustMethod = "BH", contrast=c("State","After","Before"))
summary(res_D1_t)
D1_taxo <- as.data.frame(res_D1_t) %>% rownames_to_column("MAG") %>% mutate(Drought = "D1")
#8 up, 2 down genes, with alpha = 0.05 and filtering

#DEM Second drought
counts_d2 <- counts(dds_D2) %>% data.frame() %>% rownames_to_column("Query")

df_D2 <- counts_d2 %>% left_join(eggnog_results, by = join_by(Query == query)) %>% filter(!is.na(MAGs)) %>%  
  group_by(MAGs) %>% 
  summarise(across(starts_with("RNA"), sum)) %>% 
  ungroup() %>% 
  column_to_rownames("MAGs") %>% 
  round() %>% 
  as.matrix()

D2.raw <- DESeqDataSetFromMatrix(countData = df_D2,
                                 colData = RNA_metadata_D2,
                                 design = design)
D2_obj <- DESeq(D2.raw)
res_D2_t <- results(D2_obj, alpha = 0.05, pAdjustMethod = "BH", contrast=c("State","After","Before"))
summary(res_D2_t)
D2_taxo <- as.data.frame(res_D2_t) %>% rownames_to_column("MAG") %>% mutate(Drought = "D2")
#11 up, 12 down genes, with alpha = 0.05 and filtering

#DEM Third drought
counts_d3 <- counts(dds_D3) %>% data.frame() %>% rownames_to_column("Query")

df_D3 <- counts_d3 %>% left_join(eggnog_results, by = join_by(Query == query)) %>% filter(!is.na(MAGs)) %>%  
  group_by(MAGs) %>% 
  summarise(across(starts_with("RNA"), sum)) %>% 
  ungroup() %>% 
  column_to_rownames("MAGs") %>% 
  round() %>% 
  as.matrix()

D3.raw <- DESeqDataSetFromMatrix(countData = df_D3,
                                 colData = RNA_metadata_D3,
                                 design = design)
D3_obj <- DESeq(D3.raw)
res_D3_t <- results(D3_obj, alpha = 0.05, pAdjustMethod = "BH", contrast=c("State","After","Before"))
summary(res_D3_t)
D3_taxo <- as.data.frame(res_D3_t) %>% rownames_to_column("MAG") %>% mutate(Drought = "D3")
#1 up, 0 down genes, with alpha = 0.05 and filtering

DEMAGs_all <- rbind(D1_taxo, D2_taxo, D3_taxo) %>%
  left_join(MAGs_taxo, by = join_by(MAG == MAGs)) %>% 
  mutate(DA = case_when((log2FoldChange >= 2 & padj <= 0.05) ~ "Up-regulated",
                        (log2FoldChange <= -2 & padj <= 0.05) ~ "Down-regulated",
                        TRUE ~ "No change")) %>%
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum))
  
DEMAGs_sign <-DEMAGs_all %>% 
  filter(padj <= 0.05) %>% 
  filter(log2FoldChange < -2 | log2FoldChange > 2)
write.csv(DEMAGs_sign, "Results/DEMAGs_sign.csv", row.names = FALSE)


##Visualization of DEG
Drought <- rbind(D1, D2, D3)

Drought_df <- Drought %>% 
  mutate(DA = case_when((log2FoldChange >= 2 & padj <= 0.05) ~ "Up-regulated",
                        (log2FoldChange <= -2 & padj <= 0.05) ~ "Down-regulated",
                        TRUE ~ "No change")) %>% 
  left_join(eggnog_results, by = join_by(Gene == query)) %>% 
  left_join(MAGs_taxo) %>% 
  select(Gene, log2FoldChange, padj, Drought, DA, MAGs, KEGG_ko, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(!is.na(MAGs)) %>% 
  mutate(PhylumColor = ifelse(DA == "No change", "No", as.character(Phylum)),
         DA_Period = ifelse(DA == "No change", "No change", paste(DA, Drought)),
         Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum))
write.csv(Drought_df, "Results/Drought_df.csv", row.names = FALSE)

Drought_df_sign <- Drought_df %>% 
  filter(DA != "No change")
write.csv(Drought_df_sign, "Results/Drought_df_sign.csv", row.names = FALSE)

dup_DA_patterns <- Drought_df_sign %>%
  group_by(Gene) %>%
  filter(n() > 1) %>%
  summarise(DA_pattern = paste(sort(unique(DA)), collapse = " + "), .groups = "drop") %>%
  group_by(DA_pattern) %>%
  summarise(n_genes = n())
dup_DA_patterns

Drought_df$Phylum <- factor(Drought_df$Phylum, levels = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", 
                                                          "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", 
                                                          "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota"))
Drought_df_sign$Phylum<- factor(Drought_df_sign$Phylum, levels = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", 
                                                              "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", 
                                                              "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota"))
Drought_df$DA <- factor(Drought_df$DA, levels = c("No change", "Up-regulated", "Down-regulated"))

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

Label_droughts = as_labeller(c("D1" = "Drought 1 (6h)", "D2" = "Drought 2 (24h)", "D3" = "Drought 3 (24h)"))

volcano_taxo <- Drought_df %>%
  arrange(!(Phylum %in% c("Cyanobacteriota", "Stramenopiles"))) %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  theme_classic() +
  geom_point(
    data = ~ subset(., DA == "No change"),
    aes(alpha = DA),
    fill = "lightgray",
    size = 3,
    shape = 21,
    stroke = 0.1,
    show.legend = FALSE) +
  geom_jitter(
    data = ~ subset(., DA != "No change"),
    aes(fill = Phylum, alpha = DA),
    size = 4,
    stroke = 0.1,
    shape = 21) +
  geom_col(data = Drought_df_sign,
           aes(x = -Inf, y = -Inf, fill = Phylum),
           inherit.aes = FALSE,
           col = "black",
           show.legend = TRUE) +
  scale_color_manual(values = ColorPhyla, guide = "none") +
  scale_fill_manual(values = ColorPhyla) +
  scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "No change" = 0.1), guide = "none") +
  facet_wrap(~Drought, scales = "free_x", labeller = Label_droughts) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 15), linetype="dashed", color = "lightgray", linewidth = 0.8) +
  labs(x = expression(bold("log" [2] ~ "Fold-Change")), y = expression(bold(-log ~ italic(p)*"-value"))) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.justification = "left",
        legend.spacing = unit(0, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0.2))) +
  xlim(-12,12) +
  ylim(0, 15)
volcano_taxo

Figure3 <- volcano_taxo
Figure3
ggsave("Figures/Figure3.pdf", Figure3, width = 18.2, height = 9.1)


#Info taxonomy
DEG_taxo_info <- Drought_df_sign %>%
  add_count(DA_Period, Phylum, name = "sumDEGphyla") %>%
  group_by(DA_Period) %>%
  distinct(DA_Period, Phylum, sumDEGphyla) %>% 
  mutate(sumDEG = sum(sumDEGphyla), rel_abundance = (sumDEGphyla / sumDEG)*100) %>%
  ungroup()

##Visualization of DEM
DEMAGs_sign$Phylum <- factor(DEMAGs_sign$Phylum, levels = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", 
                                                            "Patescibacteria", "Chloroflexota", "Spirochaetota", "Obazoa"))

volcano_taxo_activity <- DEMAGs_all %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  theme_classic() +
  geom_point(
    data = ~ subset(., DA == "No change"),
    fill = "lightgray",
    size = 4,
    shape = 21,
    stroke = 0.1,
    alpha = 0.5,
    show.legend = FALSE) +
  geom_point(
    data = ~ subset(., DA != "No change"),
    aes(fill = Phylum),
    size = 4,
    stroke = 0.1,
    shape = 21) +
  geom_col(data = DEMAGs_sign,
           aes(x = -Inf, y = -Inf, fill = Phylum),
           inherit.aes = FALSE,
           col = "black",
           show.legend = TRUE) +
  geom_text_repel(
    data = ~ subset(., DA != "No change"),
    aes(label = Family),
    size = 4,
    fontface = "italic",
    box.padding = 0.4,
    show.legend = FALSE,
    point.padding = 0.3) +
  scale_color_manual(values = ColorPhyla, guide = "none") +
  scale_fill_manual(values = ColorPhyla) +
  scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "No change" = 0.1), guide = "none") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 15), linetype="dashed", color = "lightgray", linewidth = 0.8) +
  facet_wrap(~Drought, scales = "free_x", labeller = Label_droughts) +
  labs(x = expression(bold("log" [2] ~ "Fold-Change")), y = expression(bold(-log ~ italic(p)*"-value"))) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 16, face="bold"),
        legend.justification = "left",
        legend.spacing = unit(0, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0.2))) +
  xlim(-5,5) +
  ylim(0, 15)
volcano_taxo_activity

#volcano_taxo_virus from DESEQ2_virus script.

Figure4 <- ggarrange(volcano_taxo_activity, volcano_taxo_virus, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), heights = c(1, 0.7), align = "v") 
Figure4
ggsave("Figures/Figure4.pdf", Figure4, width = 18.2, height = 12)
