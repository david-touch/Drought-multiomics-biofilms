

###MAG-specific scaling differential expression analysis

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
Query <- read.delim("Data/Query.tsv", sep="\t")
eggnog_results_bac <- read.delim("Data/all_annotations_bac.tsv", sep="\t")
eggnog_results_euk <- read.delim("Data/all_annotations_euk.tsv", sep="\t")
eggnog_results <- rbind(eggnog_results_bac, eggnog_results_euk)
rm(eggnog_results_bac, eggnog_results_euk)
MAGs_taxo <- read.csv("Data/MAGs_taxo_HQeuk.csv") %>% 
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum))

#Metadata preparation
RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample)
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought)
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

samples_D1 <- paste0("RNA", 1:6)
samples_D2 <- paste0("RNA", 7:12)
samples_D3 <- paste0("RNA", 13:18)

files <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata$Sample, "abundance.h5")
names(files) <- paste0("RNA", 1:18)
txi.kallisto_all <- tximport(files, type = "kallisto", txOut = TRUE)

counts_all <- txi.kallisto_all$counts
counts_all <- data.frame(counts_all) 
counts_all$Query <- rownames(counts_all)
rm(txi.kallisto_all)

##Here, load "Function" available at the end of the script.

##First drought
counts_D1 <- counts_all[,samples_D1]
filtered_D1 <- counts_D1[rowSums(counts_D1) >= 10 & rowSums(counts_D1 > 0) >= 2,]

counts_D1_taxon <- filtered_D1 %>% 
  rownames_to_column("Query") %>% 
  left_join(Query) 

feature_taxon_D1 <- counts_D1_taxon$Genome
taxa_D1 <- unique(feature_taxon_D1)

list_taxon_matrices_D1 <- lapply(taxa_D1, function(taxon) {
  as.matrix(filtered_D1)[feature_taxon_D1 == taxon, , drop = FALSE]
})
names(list_taxon_matrices_D1) <- taxa_D1

list_taxon_matrices_D1 <- Filter(Negate(is.null), list_taxon_matrices_D1)
names(list_taxon_matrices_D1) <- taxa_D1

cond_vec_D1 <- c(rep('B1', 3), rep('A1', 3)) 

results_list_D1 <- lapply(list_taxon_matrices_D1, function(taxon_mat) {
  DESeq2.result(taxon_mat, cond_vec_D1)
})
names(results_list_D1) <- taxa_D1

DEG_taxon_scale_D1 <- bind_rows(
  lapply(names(results_list_D1), function(genome) {
    res <- results_list_D1[[genome]]
    if (!is.null(res)) {
      as.data.frame(res) %>%
        mutate(Drought = "D1", MAGs = genome)
    } else {
      NULL}}))

##Second drought
counts_D2 <- counts_all[,samples_D2]
filtered_D2 <- counts_D2[rowSums(counts_D2) >= 10 & rowSums(counts_D2 > 0) >= 2,]

counts_D2_taxon <- filtered_D2 %>% 
  rownames_to_column("Query") %>% 
  left_join(Query) 

feature_taxon_D2 <- counts_D2_taxon$Genome
taxa_D2 <- unique(feature_taxon_D2)

list_taxon_matrices_D2 <- lapply(taxa_D2, function(taxon) {
  as.matrix(filtered_D2)[feature_taxon_D2 == taxon, , drop = FALSE]
})
names(list_taxon_matrices_D2) <- taxa_D2

cond_vec_D2 <- c(rep('B2', 3), rep('A2', 3)) 
 
results_list_D2 <- lapply(list_taxon_matrices_D2, function(taxon_mat) {
  DESeq2.result(taxon_mat, cond_vec_D2)
})
names(results_list_D2) <- taxa_D2

DEG_taxon_scale_D2 <- bind_rows(
  lapply(names(results_list_D2), function(genome) {
    res <- results_list_D2[[genome]]
    if (!is.null(res)) {
      as.data.frame(res) %>%
        mutate(Drought = "D2", MAGs = genome)
    } else {
      NULL}}))

##Third drought
counts_D3 <- counts_all[,samples_D3]
filtered_D3 <- counts_D3[rowSums(counts_D3) >= 10 & rowSums(counts_D2 > 0) >= 2,]

counts_D3_taxon <- filtered_D3 %>% 
  rownames_to_column("Query") %>% 
  left_join(Query) 

feature_taxon_D3 <- counts_D3_taxon$Genome
taxa_D3 <- unique(feature_taxon_D3)

list_taxon_matrices_D3 <- lapply(taxa_D3, function(taxon) {
  as.matrix(filtered_D3)[feature_taxon_D3 == taxon, , drop = FALSE]
})
names(list_taxon_matrices_D3) <- taxa_D3

cond_vec_D3 <- c(rep('B3', 3), rep('A3', 3)) 

results_list_D3 <- lapply(list_taxon_matrices_D3, function(taxon_mat) {
  DESeq2.result(taxon_mat, cond_vec_D3)
})
names(results_list_D3) <- taxa_D3

DEG_taxon_scale_D3 <- bind_rows(
  lapply(names(results_list_D3), function(genome) {
    res <- results_list_D3[[genome]]
    if (!is.null(res)) {
      as.data.frame(res) %>%
        mutate(Drought = "D3", MAGs = genome)
    } else {
      NULL}}))


##Combine DEG taxon-scaled
DEG_taxon_scaling <- rbind(DEG_taxon_scale_D1, DEG_taxon_scale_D2, DEG_taxon_scale_D3) %>% 
  rownames_to_column("Query") %>% 
  mutate(DA = case_when((log2FoldChange >= 2 & padj <= 0.05) ~ "Up-regulated",
                        (log2FoldChange <= -2 & padj <= 0.05) ~ "Down-regulated",
                        TRUE ~ "No change")) %>% 
  left_join(select(eggnog_results, -MAGs), by = join_by(Query == query)) %>% 
  left_join(MAGs_taxo) %>% 
  select(Query, log2FoldChange, padj, Drought, DA, MAGs, KEGG_ko, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(!is.na(MAGs)) %>% 
  mutate(PhylumColor = ifelse(DA == "No change", "No", as.character(Phylum)),
         DA_Period = ifelse(DA == "No change", "No change", paste(DA, Drought)))

DEG_taxon_scaling_sign <- DEG_taxon_scaling %>%
  filter(DA != "No change")
write.csv(DEG_taxon_scaling_sign, "Results/DEG_taxon_scaling_sign.csv", row.names = FALSE)

dup_DA_patterns <- DEG_taxon_scaling_sign %>%
  group_by(Query) %>%
  filter(n() > 1) %>%
  summarise(DA_pattern = paste(sort(unique(DA)), collapse = " + "), .groups = "drop") %>%
  group_by(DA_pattern) %>%
  summarise(n_genes = n())
dup_DA_patterns

DEG_taxon_scaling_sign$Phylum<- factor(DEG_taxon_scaling_sign$Phylum, levels = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", 
                                                                                 "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", 
                                                                                 "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota"))

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


##Plotting
volcano_taxo_taxon_scale <- DEG_taxon_scaling %>%
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
  geom_col(data = DEG_taxon_scaling_sign,
           aes(x = -Inf, y = -Inf, fill = Phylum),
           inherit.aes = FALSE,
           col = "black",
           show.legend = TRUE) +
  scale_color_manual(values = ColorPhyla, guide = "none") +
  scale_fill_manual(values = ColorPhyla) +
  scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "No change" = 0.1), guide = "none") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 20), linetype="dashed", color = "lightgray", linewidth = 0.8) +
  facet_wrap(~Drought, scales = "free_x", labeller = Label_droughts) +
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
  xlim(-13,13) +
  ylim(0, 20)
volcano_taxo_taxon_scale

FigureS10 <- volcano_taxo_taxon_scale
FigureS10 
ggsave("Figures/FigureS10.pdf", FigureS10, width = 18.2, height = 9.1)

#Info taxonomy
DEG_taxon_taxo_info_1 <- DEG_taxon_scaling_sign %>%
  add_count(Phylum, name = "sumDEGphyla") %>%
  distinct(Phylum, sumDEGphyla) %>% 
  mutate(sumDEG = sum(sumDEGphyla), rel_abundance = (sumDEGphyla / sumDEG)*100) %>%
  ungroup()
DEG_taxon_taxo_info_2 <- DEG_taxon_scaling_sign %>%
  add_count(Drought, name = "sumDEGdrought") %>%
  distinct(Drought, sumDEGdrought) %>% 
  mutate(sumDEG = sum(sumDEGdrought), rel_abundance = (sumDEGdrought / sumDEG)*100) %>%
  ungroup()


##Function
DESeq2.result <- function(Xmat, cond) {
  Xmat = round(Xmat)
  storage.mode(Xmat) <- 'integer'
  colData <- data.frame(condition = cond)
  
  # Filter out rows with all zero counts
  Xmat <- Xmat[rowSums(Xmat) > 0, , drop = FALSE]
  if (nrow(Xmat) < 2) return(NULL)
  
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)
  colData(dds)$condition <- factor(colData(dds)$condition, levels = unique(cond))
  
  tryCatch({
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    res <- results(dds)
    return(res)
  }, error = function(e) {
    message("Fallback for dispersion fitting on current taxon: using gene-wise estimates with poscounts")
    
    # Safe fallback approach
    tryCatch({
      dds <- estimateSizeFactors(dds, type = "poscounts")
      
      dds <- estimateDispersionsGeneEst(dds)
      dispGeneEst <- mcols(dds)$dispGeneEst
      
      if (is.null(dispGeneEst)) {
        message("  -> dispGeneEst is NULL, skipping.")
        return(NULL)
      }
      
      # Keep only rows with valid dispersion estimates
      valid_idx <- which(!is.na(dispGeneEst))
      if (length(valid_idx) < 2) {
        message("  -> Too few valid dispersion estimates (", length(valid_idx), "). Skipping.")
        return(NULL)
      }
      
      # Subset both dds and dispersion estimates
      dds <- dds[valid_idx, ]
      dispGeneEst <- dispGeneEst[valid_idx]
      
      # Final safety check
      if (length(dispGeneEst) != nrow(dds)) {
        message("  -> Mismatch after filtering: ", length(dispGeneEst), " vs ", nrow(dds))
        return(NULL)
      }
      
      if (any(is.na(dispGeneEst))) {
        message("  -> Still NA values after filtering! Skipping.")
        return(NULL)
      }
      
      # Assign dispersions and run Wald test
      names(dispGeneEst) <- rownames(dds)
      dispersions(dds) <- dispGeneEst
      
      dds <- nbinomWaldTest(dds)
      res <- results(dds)
      message("  -> Successfully computed DESeq2 results for this taxon")
      return(res)
      
    }, error = function(e2) {
      message("  -> Fallback failed: ", conditionMessage(e2))
      return(NULL)
    })
  })
}
