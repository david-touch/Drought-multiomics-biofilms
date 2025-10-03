

###Visualisation of differentially abundant MAGs - taxon-scaled

#Load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(clusterProfiler)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data
Drought_TS_gene_sign <- read.csv("Results/DEG_taxon_scaling_sign.csv")

DEGs_TS_func_taxo <- Drought_TS_gene_sign %>% 
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:"),
         KEGG_ko = na_if(KEGG_ko, "-"))
  

#Per drought and per group df
D1_df <- DEGs_TS_func_taxo %>% 
  filter(Drought == "D1")

D1_df_sign_up_cyano <- D1_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Cyanobacteriota")
D1_df_sign_down_cyano <- D1_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Cyanobacteriota")
D1_df_sign_up_euk <- D1_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Stramenopiles")
D1_df_sign_down_euk <- D1_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Stramenopiles")
D1_df_sign_up_hetero <- D1_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")
D1_df_sign_down_hetero <- D1_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")

D2_df <- DEGs_TS_func_taxo %>% 
  filter(Drought == "D2")

D2_df_sign_up_cyano <- D2_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Cyanobacteriota")
D2_df_sign_down_cyano <- D2_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Cyanobacteriota")
D2_df_sign_up_euk <- D2_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Stramenopiles")
D2_df_sign_down_euk <- D2_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Stramenopiles")
D2_df_sign_up_hetero <- D2_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")
D2_df_sign_down_hetero <- D2_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")


D3_df <- DEGs_TS_func_taxo %>% 
  filter(Drought == "D3")

D3_df_sign_up_cyano <- D3_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Cyanobacteriota")
D3_df_sign_down_cyano <- D3_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Cyanobacteriota")
D3_df_sign_up_euk <- D3_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum == "Stramenopiles")
D3_df_sign_down_euk <- D3_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum == "Stramenopiles")
D3_df_sign_up_hetero <- D3_df %>% filter(log2FoldChange >= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")
D3_df_sign_down_hetero <- D3_df %>% filter(log2FoldChange <= 2) %>% filter(Phylum != "Cyanobacteriota" & Kingdom != "Eukaryote")


##Enrichment analysis

#Cyanobacteria
D1_geneUp_cyano <- D1_df_sign_up_cyano$KEGG_ko
D1_geneDown_cyano <- D1_df_sign_down_cyano$KEGG_ko
D1_kk_up_cyano <- enrichKEGG(gene = D1_geneUp_cyano,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.01)
D1_kk_up_df_cyano <- data.frame(D1_kk_up_cyano) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Up")
#No change up D1 cyano
D1_kk_down_cyano <- enrichKEGG(gene = D1_geneDown_cyano,
                               organism     = 'ko',
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.01)
#No change down D1 cyano

D2_geneUp_cyano <- D2_df_sign_up_cyano$KEGG_ko
D2_geneDown_cyano <- D2_df_sign_down_cyano$KEGG_ko
D2_kk_up_cyano <- enrichKEGG(gene = D2_geneUp_cyano,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
D2_kk_up_df_cyano <- data.frame(D2_kk_up_cyano) %>% filter(category != "Human Diseases", category != "Organismal Systems", Count > 5) %>% mutate(DA = "Up", Taxonomy = "Cyanobacteriota", Drought = "D2")
#Enrichment up D2 cyano
D2_kk_down_cyano <- enrichKEGG(gene = D2_geneDown_cyano,
                               organism     = 'ko',
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
D2_kk_down_df_cyano <- data.frame(D2_kk_down_cyano) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Down", Taxonomy = "Cyanobacteriota", Drought = "D2")
#No change down D2 cyano

D3_geneUp_cyano <- D3_df_sign_up_cyano$KEGG_ko
D3_geneDown_cyano <- D3_df_sign_down_cyano$KEGG_ko
D3_kk_up_cyano <- enrichKEGG(gene = D3_geneUp_cyano,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
#No change up D3 cyano
D3_kk_down_cyano <- enrichKEGG(gene = D3_geneDown_cyano,
                               organism     = 'ko',
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
#No change down D3 cyano


#Stramenopiles
D1_geneUp_euk <- D1_df_sign_up_euk$KEGG_ko
D1_geneDown_euk <- D1_df_sign_down_euk$KEGG_ko
D1_kk_up_euk <- enrichKEGG(gene = D1_geneUp_euk,
                           organism     = 'ko',
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01)
D1_kk_up_df_euk <- data.frame(D1_kk_up_euk) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Up")
#No change up D1 euk
D1_kk_down_euk <- enrichKEGG(gene = D1_geneDown_euk,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.01)
D1_kk_down_df_euk <- data.frame(D1_kk_down_euk) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Down", Taxonomy = "Stramenopiles", Drought = "D1")
#Enrichement down D1 euk

D2_geneUp_euk <- D2_df_sign_up_euk$KEGG_ko
D2_geneDown_euk <- D2_df_sign_down_euk$KEGG_ko
D2_kk_up_euk <- enrichKEGG(gene = D2_geneUp_euk,
                           organism     = 'ko',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
D2_kk_up_df_euk <- data.frame(D2_kk_up_euk) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Up", Taxonomy = "Stramenopiles", Drought = "D2")
#Enrichement up D2 euk
D2_kk_down_euk <- enrichKEGG(gene = D2_geneDown_euk,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
D2_kk_down_df_euk <- data.frame(D2_kk_down_euk) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Down", Taxonomy = "Stramenopiles", Drought = "D2")
#No change up D2 euk

D3_geneUp_euk <- D3_df_sign_up_euk$KEGG_ko
D3_geneDown_euk <- D3_df_sign_down_euk$KEGG_ko
D3_kk_up_euk <- enrichKEGG(gene = D3_geneUp_euk,
                           organism     = 'ko',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
#No change up D3 euk
D3_kk_down_euk <- enrichKEGG(gene = D3_geneDown_euk,
                             organism     = 'ko',
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
#No change down D3 euk


#Heterotroph
D1_geneUp_hetero <- D1_df_sign_up_hetero$KEGG_ko
D1_geneDown_hetero <- D1_df_sign_down_hetero$KEGG_ko
D1_kk_up_hetero <- enrichKEGG(gene = D1_geneUp_hetero,
                              organism     = 'ko',
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.01)
D1_kk_up_df_hetero <- data.frame(D1_kk_up_hetero) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Up", Taxonomy = "Heterotroph", Drought = "D1")
#No change up D1
D1_kk_down_hetero <- enrichKEGG(gene = D1_geneDown_hetero,
                                organism     = 'ko',
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.01)
D1_kk_down_df_hetero <- data.frame(D1_kk_down_hetero) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Down")
#No change down D1 hetero

D2_geneUp_hetero <- D2_df_sign_up_hetero$KEGG_ko
D2_geneDown_hetero <- D2_df_sign_down_hetero$KEGG_ko
D2_kk_up_hetero <- enrichKEGG(gene = D2_geneUp_hetero,
                              organism     = 'ko',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)
D2_kk_up_df_hetero <- data.frame(D2_kk_up_hetero) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Up", Taxonomy = "Heterotroph", Drought = "D2")
#No change up D2 hetero
D2_kk_down_hetero <- enrichKEGG(gene = D2_geneDown_hetero,
                                organism     = 'ko',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05)
D2_kk_down_df_hetero <- data.frame(D2_kk_down_hetero) %>% filter(category != "Human Diseases", Count > 5) %>% mutate(DA = "Down")
#No change down D2 hetero

D3_geneUp_hetero <- D3_df_sign_up_hetero$KEGG_ko
D3_geneDown_hetero <- D3_df_sign_down_hetero$KEGG_ko
D3_kk_up_hetero <- enrichKEGG(gene = D3_geneUp_hetero,
                              organism     = 'ko',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05)
#No change up D3 hetero
D3_kk_down_hetero <- enrichKEGG(gene = D3_geneDown_hetero,
                                organism     = 'ko',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05)
#No change down D3 hetero

Enrichment_results_TS <- rbind(D2_kk_up_df_cyano, D1_kk_down_df_euk, D2_kk_up_df_euk, D2_kk_down_df_euk, D1_kk_up_df_hetero, D2_kk_up_df_hetero)
write.csv(Enrichment_results_TS, "Results/Enrichment_results_TS.csv")
