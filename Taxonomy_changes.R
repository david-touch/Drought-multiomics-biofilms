

###Taxonomic changes on MAGs - metagenomic and metatranscpritomic 

#Load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)
library(purrr)
library(broom)
library(ggh4x)
library(rstatix)
library(ggtext)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)

#Load data
DNA_metadata <- read_csv("Data/DNA_metadata.csv") %>% data.frame()
Rel_abundance_DNA_samples <- read.delim("Data/DNA_coverage_samples.tsv", sep="\t")
MAGs_taxo <- read.csv("Data/MAGs_taxo_HQeuk.csv")

#Create a phyloseq object
otu_table_DNA <- Rel_abundance_DNA_samples %>% column_to_rownames("Genome") %>% otu_table(taxa_are_rows = T)
sample_data_DNA <- DNA_metadata %>% column_to_rownames("Name") %>% mutate(ID = Sample) %>% sample_data()
tax_table_DNA <- MAGs_taxo %>% filter(MAGs %in% Rel_abundance_DNA_samples$Genome) %>% column_to_rownames("MAGs") %>% 
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum)) %>% as.matrix() %>%  tax_table()

DNA_phylo <- merge_phyloseq(otu_table_DNA, sample_data_DNA, tax_table_DNA)

DNA_phylo_rel  <- transform_sample_counts(DNA_phylo, function(x) x / sum(x))

#Plotting
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

Phyla_order <- c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Stramenopiles", "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota")

Label_droughts = as_labeller(c("1" = "Drought 1 (6h)", "2" = "Drought 2 (24h)", "3" = "Drought 3 (24h)"))

Label_samples_DNA <- c("S4" = "Before", "S5" = "After",
                       "S8" = "Before", "S9" = "After",
                       "S17" = "Before", "S18" = "After")

DNA_barplot_droughts <- DNA_phylo_rel %>%
  ps_filter(!is.na(Period)) %>% 
  ps_arrange(Sample) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "ID",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9,
    palette = ColorPhyla) + 
  theme_classic() +
  xlab("")+
  ylab("Relative abundance\n(gDNA)") +
  scale_x_discrete(labels = Label_samples_DNA) +
  facet_wrap(~Period, scales = "free_x", labeller = Label_droughts) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"))
DNA_barplot_droughts

DNA_barplot_all <- DNA_phylo_rel %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9,
    palette = ColorPhyla) + 
  theme_classic() +
  xlab("Days of growth")+
  ylab("Relative abundance\n(gDNA)") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"))
DNA_barplot_all

#NMDS taxomony
DNA_phylo_hell <- microbiome::transform(DNA_phylo, 'hell')
matrix_dist <- as.matrix(t(data.frame(otu_table(DNA_phylo_hell))))
metadata_dist <- data.frame(sample_data(DNA_phylo_hell)) %>% rownames_to_column("Samples")
dist_MAGs <- vegdist(matrix_dist, method = "bray")
nmds_MAGs <- metaMDS(dist_MAGs, trymax=100)
stress <- nmds_MAGs$stress
stress

scores_nmds_MAGs <- scores(nmds_MAGs) %>%
  data.frame() %>% 
  rownames_to_column("Name") %>% 
  left_join(DNA_metadata)

env <- data.frame(sample_data(DNA_phylo_hell)) %>% 
  rownames_to_column("Samples") %>% 
  select(-c("Samples", "Sample", "Time", "Period", "ChlA1", "ChlA2", "ChlA3", "BA1", "BA2", "BA3", "ID"))
set.seed(19950930)
fit <- envfit(nmds_MAGs, env, permutations = 999, na.rm = TRUE)

fit_vec <- data.frame(fit$vectors$arrows * sqrt(fit$vectors$r), P = fit$vectors$pvals)
fit_vec$Environment <- rownames(fit_vec)
fit_fac <- data.frame(fit$factors$centroids * sqrt(fit$factors$r), P = fit$factors$pvals)
fit_fac$Environment <- rownames(fit_fac) 
fit_fac <- fit_fac %>% slice(2L) %>% mutate(Environment = "Drought")
fit_env <- rbind(fit_vec,fit_fac)
fit_env_pv <- subset(fit_env, P<=0.05) %>%
  group_by(Environment) %>%
  mutate(EnvironmentSign = case_when(
    Environment %in% c("Succession", "Drought") ~ paste0("\"", Environment, "*\""),
    TRUE ~ paste0("\"", Environment, "\"")
  )) %>%
  mutate(EnvironmentSign = as.character(EnvironmentSign))

NMDS_MAGs <- scores_nmds_MAGs %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(color = "darkgray", size = 3) +
  theme_classic() +
  geom_segment(data=scores_nmds_MAGs, aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.5, color = "#b03a2e") +
  labs(caption="Stress 0.03") +
  theme(axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=1, vjust = 22),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent")) +
  geom_segment(data = fit_env_pv, aes(x = 0, xend = NMDS1 * 0.4, y = 0, yend = NMDS2 * 0.4),
               arrow = arrow(length = unit(0.3, "cm")), colour = "darkgrey", lwd = 0.75) +
  ggrepel::geom_text_repel(data = fit_env_pv,
                           aes(x = NMDS1 * 0.4, y = NMDS2 * 0.4,
                               label = EnvironmentSign),
                           parse = TRUE,
                           cex = 5, direction = "x", colour = "black",
                           segment.size = 0.25,
                           nudge_x = c(0.01, 0.01), nudge_y = c(0.01, 0.01))
NMDS_MAGs

FigureS4 <- ggarrange(DNA_barplot_all, NMDS_MAGs, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), heights = c(1, 0.7), legend="top")
FigureS4
ggsave("Figures/FigureS4.pdf", FigureS4, width = 18.2, height = 18.2)

adonis_bc <- adonis2(matrix_dist ~ Drought+Succession, data=metadata_dist, permutations=999, method="bray", by = "terms")
adonis_bc


##Metatranscriptomic taxonomy

#load data
Query <- read.delim("Data/Query.tsv", sep="\t")

RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample)
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought, levels = c("B1","A1", "B2", "A2", "B3", "A3"))
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

#All droughts
files <- file.path(".", "Data", "Kallisto", "HQ_Genes", RNA_metadata$Sample, "abundance.tsv")
names(files) <- paste0("RNA", 1:18)

tpm_list <- imap(files, ~ {
  read_tsv(.x, show_col_types = FALSE) %>%
    select(target_id, !!.y := tpm)  # Rename 'tpm' column to sample name
})
tpm_combined <- purrr::reduce(tpm_list, full_join, by = "target_id")

tpm_filtered <- tpm_combined[rowSums(tpm_combined >= 0.1) >= 2, ]


##Create a phyloseq object
otu_table_RNA <- tpm_filtered %>% column_to_rownames("target_id") %>% otu_table(taxa_are_rows = T)
sample_data_RNA <- RNA_metadata %>% select(-Name) %>% mutate(ID = Sample) %>% column_to_rownames("Sample") %>% sample_data()
tax_table_RNA <- Query %>% left_join(MAGs_taxo, by = join_by(Genome == MAGs)) %>% column_to_rownames("Query") %>% 
  mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum)) %>% 
  select(-Genome) %>% as.matrix() %>% tax_table()

RNA_phylo <- merge_phyloseq(otu_table_RNA, sample_data_RNA, tax_table_RNA)

RNA_phylo_rel  <- transform_sample_counts(RNA_phylo, function(x) x / sum(x))

Phyla <- tax_glom(RNA_phylo_rel, taxrank = "Phylum") %>% 
  ps_melt() %>% 
  group_by(Phylum) %>% 
  summarise(av_abund = mean(Abundance)) %>%
  arrange(av_abund)

Active_phyla <- Phyla %>% pull(Phylum)
Top_active_phyla <- Phyla %>% filter(av_abund >= 0.02) %>% pull(Phylum)

Label_samples_RNA <- c("B1" = "Before", "A1" = "After",
                   "B2" = "Before", "A2" = "After",
                   "B3" = "Before", "A3" = "After")

Label_droughts_RNAall = as_labeller(c("1" = "Drought 1 (6h)", "2" = "Drought 2 (24h)", "3" = "Drought 3 (24h)"))

#Plotting
RNA_barplot_droughts <- RNA_phylo_rel %>%
  ps_filter(!is.na(Period)) %>% 
  merge_samples("Drought", fun = "mean") %>% 
  ps_arrange(State) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "SAMPLE",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9, 
    palette = ColorPhyla) + 
  theme_classic() +
  ylab("Relative abundance\n(mRNA)") +
  scale_x_discrete(labels = Label_samples_RNA) +
  facet_wrap(~Period, scales = "free_x", labeller = Label_droughts) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_blank())
RNA_barplot_droughts

#Relative abundance
Kingdom_RNA <- tax_glom(RNA_phylo_rel, taxrank = "Kingdom") %>% 
  ps_melt() %>% 
  group_by(Kingdom) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
Phylum_RNA <- tax_glom(RNA_phylo_rel, taxrank = "Phylum") %>% 
  ps_melt() %>% 
  group_by(Phylum) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))


#TPM differences
df_transcripts <- ps_melt(RNA_phylo) %>%
  group_by(Sample, Phylum) %>% 
  mutate(SumPhyla = sum(Abundance),
         Phyla_drought = paste(Phylum, Drought)) %>% 
  distinct(Sample, Phylum, .keep_all = T) %>% 
  ungroup() %>% 
  mutate(Phylum = factor(Phylum, levels = Active_phyla))
shapiro.test(df_transcripts$SumPhyla) #Not normal

df_transcripts$State <- factor(df_transcripts$State, levels = c("After", "Before"), labels = c("Post-Drought", "Pre-Drought"))
df_transcripts$Drought <- factor(df_transcripts$Drought, levels = c("B1", "A1", "B2", "A2", "B3", "A3"))

effect_sizes_all <- df_transcripts %>%
  group_by(Period, Phylum) %>%
  wilcox_effsize(SumPhyla ~ Drought, paired = FALSE, 
                 effect.size = "r", 
                 alternative = "two.sided") %>%
  mutate(label = paste0("r = ", sprintf("%.2f", effsize))) %>% 
  mutate(label = ifelse(abs(effsize) >= 0.75,
                        paste0("r = **", sprintf("%.2f", effsize), "**"),  # bold for large effect
                        paste0("r = ", sprintf("%.2f", effsize))))

effect_sizes_top <- effect_sizes_all %>% filter(Phylum %in% Top_active_phyla)

Transcript_top <- df_transcripts %>%
  filter(Phylum %in% Top_active_phyla) %>% 
  ggplot(aes(y = SumPhyla, x = Phylum)) +
  geom_boxplot(aes(fill = Drought, group = State)) +
  theme_classic() +
  scale_fill_manual(values = c(
    "B1" = "#a6cee3",
    "B2" = "#1f78b4",
    "B3" = "#295184",
    "A1" = "#f1948a",
    "A2" = "#ec7063",
    "A3" = "#b03a2e"),
    breaks = c("B2", "A2"),
    labels = c("Pre-drought", "Post-drought")) +
  scale_y_log10(breaks = c(1000, 10000, 100000), labels = c("1K", "10K", "100K")) +
  scale_x_discrete(position = "top") +
  ylab(expression(bold("Transcript per million (log"[10]*")"))) +
  xlab("") +
  facet_grid(Phylum ~ Period, scales = "free_y", space = "free_y", switch = "y") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 14),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray80", size = 0.2),
        legend.position = "bottom",
        legend.justification = "right",
        legend.margin = margin(t = -10),
        legend.box.margin = margin(t = -10),
        plot.margin = margin(t = 0, r = 1, b = 0, l = 2, unit = "cm")) +
  stat_compare_means(aes(group = Phyla_drought), label = "p.signif", method = "kruskal.test", hide.ns = TRUE, size = 10, vjust = 1.5, hjust = 26) +
  coord_flip() +
  geom_richtext(
    data = effect_sizes_top,
    aes(x = Phylum, y = 1000, label = label),
    inherit.aes = FALSE,
    size = 5,
    hjust = 0.8,
    vjust = 1.5,
    fill = NA,
    label.color = NA)
Transcript_top

Barplot_droughts <- ggarrange(DNA_barplot_droughts, RNA_barplot_droughts, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), common.legend =  T, legend = "right")
Figure2 <- ggarrange(Barplot_droughts, Transcript_top, ncol = 1, nrow = 2, labels = c("", "C"), font.label = list(size = 18), heights = c(1, 0.75))
Figure2
ggsave("Figures/Figure2.pdf", Figure2, width = 18.2, height = 12)


FigureS6 <- RNA_phylo_rel %>%
  ps_filter(!is.na(Period)) %>% 
  ps_arrange(State) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "SAMPLE",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9, 
    palette = ColorPhyla) + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0.01)) +
  ylab("Relative abundance") +
  facet_nested(~ Period + State, scales = "free_x", space = "free_x", nest_line = element_line(linetype = 1), labeller = labeller(Period = Label_droughts)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 16, face="bold"))
FigureS6
ggsave("Figures/FigureS6.pdf", FigureS4, width = 18.2, height = 9.1)




FigureS8 <- df_transcripts %>%
  ggplot(aes(y = SumPhyla, x = Phylum, fill = Drought, group = State)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10(breaks = c(100, 1000, 10000, 100000), labels = c("100", "1K", "10K", "100K")) +
  scale_fill_manual(values = c(
    "B1" = "#a6cee3",
    "B2" = "#1f78b4",
    "B3" = "#295184",
    "A1" = "#f1948a",
    "A2" = "#ec7063",
    "A3" = "#b03a2e"),
    breaks = c("B2", "A2"),
    labels = c("Pre-drought", "Post-drought")) +
  ylab(expression(bold("Transcript per million (log"[10]*")"))) +
  facet_grid(Phylum ~ Period, scales = "free_y", space = "free_y", labeller = Label_droughts) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 14),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray80", size = 0.2)) +
  stat_compare_means(aes(group = Phyla_drought), label = "p.signif", method = "kruskal.test", hide.ns = TRUE, size = 10, vjust = 1.5, hjust = 24.5) +
  coord_flip() +
  geom_richtext(
    data = effect_sizes_all,
    aes(x = Phylum, y = 100, label = label),
    inherit.aes = FALSE,
    size = 4,
    fill = NA,
    hjust = 0.2,
    vjust = 1.5,
    label.color = NA)
FigureS8
ggsave("Figures/FigureS8.pdf", FigureS8, width = 18.2, height = 14)
