

###Differential expression analysis on Viruses

#Load libraries
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(ggpubr)
library(tximport)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data
taxo_virus_all <- read_csv("Data/taxo_virus_all.csv")

RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample)
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought)
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

#List negative control
file <- file.path(".", "Data", "Kallisto", "Virus", "RNANEG", "abundance.tsv")
names(file) <- paste0("RNANEG")

NEG_list <- imap(file, ~ {
  read_tsv(.x, show_col_types = FALSE) %>%
    select(target_id, !!.y := tpm)  # Rename 'tpm' column to sample name
})
NEG_list <- purrr::reduce(NEG_list, full_join, by = "target_id") %>%  filter(RNANEG > 0)


#DEViruses

#DEV First drought
RNA_metadata_D1 <- RNA_metadata %>% filter(Time == 4 | Time == 5)
files_D1 <- file.path(".", "Data", "Kallisto", "Virus", RNA_metadata_D1$Sample, "abundance.h5")
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
D1 <- as.data.frame(res_D1_ordered) %>% rownames_to_column("Virus") %>% mutate(Drought = "D1")
#138 up, 73 down genes, with alpha = 0.05 and filtering


#DEV Second drought
RNA_metadata_D2 <- RNA_metadata %>% filter(Time == 8 | Time == 9)
files_D2 <- file.path(".", "Data", "Kallisto", "Virus", RNA_metadata_D2$Sample, "abundance.h5")
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
D2 <- as.data.frame(res_D2_ordered) %>% rownames_to_column("Virus") %>% mutate(Drought = "D2")
#153 up, 147 down genes, with alpha = 0.05 and filtering 

#DEV Third drought
RNA_metadata_D3 <- RNA_metadata %>% filter(Time == 17 | Time == 18)
files_D3 <- file.path(".", "Data", "Kallisto", "Virus", RNA_metadata_D3$Sample, "abundance.h5")
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
D3 <- as.data.frame(res_D3_ordered) %>% rownames_to_column("Virus") %>% mutate(Drought = "D3")
#3 up, 1 down genes, with alpha = 0.05 and filtering 


Drought_DE_virus <- rbind(D1, D2, D3) 

Drought_DE_virus_df <- Drought_DE_virus %>% left_join(taxo_virus_all) %>% 
  mutate(DA = case_when((log2FoldChange >= 2 & padj <= 0.05) ~ "Up-regulated",
                        (log2FoldChange <= -2 & padj <= 0.05) ~ "Down-regulated",
                        TRUE ~ "No change"),
         DA_Period = ifelse(DA == "No change", "No change", paste(DA, Drought))) %>% 
  filter(Domain != "Domain")

write.csv(Drought_DE_virus, "Results/Drought_DE_virus_all.csv", row.names = F)

Drought_DE_virus_sign <- Drought_DE_virus_df %>% filter(DA != "No change")
write.csv(Drought_DE_virus_sign, "Results/Drought_DE_virus_sign.csv", row.names = F)

ColorVirus <- c(
  Caudoviricetes = "#8DD3C7",
  Chrymotiviricetes = "#FFED6F",
  Alsuviricetes = "#CCEBC5",
  "LTR retrotransposons" =  "#BEBADA",
  Megaviricetes = "#FB8072",
  Pokkesviricetes = "#80B1D3",
  Magsaviricetes = "#FDB462",
  Other = "lightgray")

Drought_DE_virus_df <- Drought_DE_virus_df %>%
  mutate(Class = fct_other(Class, keep = names(ColorVirus)[names(ColorVirus) != "Other"], other_level = "Other"))

Drought_DE_virus_sign <- Drought_DE_virus_sign %>%
  mutate(Class = fct_other(Class, keep = names(ColorVirus)[names(ColorVirus) != "Other"], other_level = "Other"))

Label_droughts = as_labeller(c("D1" = "Drought 1 (6h)", "D2" = "Drought 2 (24h)", "D3" = "Drought 3 (24h)"))

volcano_taxo_virus <- Drought_DE_virus_df %>%
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
    aes(fill = Class, alpha = DA),
    size = 4,
    width = 0.2,   # horizontal nudge
    height = 0.1, 
    stroke = 0.1,
    shape = 21) +
  geom_col(data = Drought_DE_virus_sign,
           aes(x = -Inf, y = -Inf, fill = Class),
           inherit.aes = FALSE,
           col = "black",
           show.legend = TRUE) +
  scale_color_manual(values = ColorVirus, guide = "none") +
  scale_fill_manual(values = ColorVirus) +
  scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "No change" = 0.1), guide = "none") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 12), linetype="dashed", color = "lightgray", linewidth = 0.8) +
  facet_wrap(~Drought, scales = "free_x", labeller = Label_droughts) +
  labs(x = expression(bold("log" [2] ~ "Fold-Change")), y = expression(bold(-log ~ italic(p)*"-value"))) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        legend.justification = "left",
        legend.spacing = unit(0, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0.2))) +
  xlim(-10,10) +
  ylim(0, 12)
volcano_taxo_virus

#
DEV_info_drought <- Drought_DE_virus_sign %>% 
  add_count(Drought, name = "sumDEVdrought")

DEV_info_taxo <- Drought_DE_virus_sign %>% 
  add_count(DA_Period, Phylum, name = "sumDEVclass") %>% 
  group_by(DA_Period) %>%
  distinct(DA_Period, Class, sumDEVclass, .keep_all = T)
