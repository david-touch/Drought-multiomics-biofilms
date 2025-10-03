

###Taxonomic analysis of DNA and RNA viruses 

#Load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(microViz)
library(RColorBrewer)
library(ggpubr)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)

#Load and prepare data
DNA_metadata <- read_csv("Data/DNA_metadata.csv") %>% data.frame()
RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

taxo_virus_all <- read_csv("Data/taxo_virus_all.csv") %>% 
  filter(!is.na(Domain)) %>% 
  column_to_rownames("Virus")

DNA_virus_table <- read.delim("Data/DNA_coverm_trimmed_mean_all_viruses.tsv", sep="\t") %>% column_to_rownames("Contig")
RNA_virus_table <- read.delim("Data/RNA_coverm_trimmed_mean_all_viruses.tsv", sep="\t") %>% column_to_rownames("Contig")

#Create DNA phyloseq object
otu_table_DNA_virus <- DNA_virus_table %>% setNames(str_extract(colnames(.), "^Sentinel_\\d+")) %>% otu_table(taxa_are_rows = T)
sample_data_DNA_virus <- DNA_metadata %>% mutate(Name = str_extract(Name, "^Sentinel_\\d+")) %>% column_to_rownames("Name") %>% mutate(ID = Sample) %>% sample_data()
tax_table_DNA_virus <- taxo_virus_all %>% as.matrix() %>% tax_table()

DNA_phylo_virus <- merge_phyloseq(otu_table_DNA_virus, sample_data_DNA_virus, tax_table_DNA_virus) %>% tax_fix(
  min_length = 4, sep = " ", anon_unique = TRUE,
  suffix_rank = "classified")

DNA_phylo_virus_rel  <- transform_sample_counts(DNA_phylo_virus, function(x) x / sum(x))
#2242 DNA vMAGs

DNA_colors <- c(
  Caudoviricetes = "#8DD3C7",
  Arfiviricetes = "#FFFFB3",
  "LTR retrotransposons" =  "#BEBADA",
  Megaviricetes = "#FB8072", 
  Pokkesviricetes = "#80B1D3",
  Other = "lightgray")

#Plotting
DNA_barplot_virus <- DNA_phylo_virus_rel %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Class", n_taxa = 5,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9, 
    palette = DNA_colors) + 
  theme_classic() +
  xlab("")+
  ylab("Relative abundance\n(DNA vMAGs)") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"))
DNA_barplot_virus

##Normalize per sequencing depth
DNA_otu <- as(otu_table(DNA_phylo_virus), "matrix")
if (taxa_are_rows(DNA_phylo_virus)) {
  DNA_otu <- t(DNA_otu)
}

DNA_depths <- sample_data(DNA_phylo_virus)$Depth
names(DNA_depths) <- rownames(sample_data(DNA_phylo_virus))
DNA_depths <- DNA_depths[rownames(DNA_otu)]
DNA_otu_norm <- sweep(DNA_otu, 1, DNA_depths, "/") * 1e6

DNA_phylo_virus_norm <- phyloseq(otu_table(DNA_otu_norm, taxa_are_rows = FALSE), tax_table(DNA_phylo_virus), sample_data(DNA_phylo_virus))


DNA_norm_class <- tax_glom(DNA_phylo_virus_norm, taxrank = "Class")
DNA_plot_df <- psmelt(DNA_norm_class)

DNA_top_classes <- DNA_plot_df %>%
  group_by(Class) %>%
  summarise(mean_ab = mean(Abundance), .groups = "drop") %>%
  slice_max(mean_ab, n = 5) %>%
  pull(Class)

DNA_plot_df <- DNA_plot_df %>%
  mutate(Class = ifelse(Class %in% DNA_top_classes, Class, "Other"))

DNA_plot_df$Class <- factor(DNA_plot_df$Class, levels = c("Caudoviricetes", "Arfiviricetes", "LTR retrotransposons", "Megaviricetes", "Pokkesviricetes", "Other"))

DNA_plot_df$Succession <- factor(DNA_plot_df$Succession, levels = c("5","10","14","19","24","26","28","33",
                                                "38","40","42","47","54","61","67","75",
                                                "82","87","89","91","96","103"))

DNA_barplot_virus_norm <- ggplot(DNA_plot_df, aes(x = Succession, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", colour = "grey5", size = 0.1, width = 0.9, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = DNA_colors) +
  theme_classic() +
  xlab("Days of growth") +
  ylab("Counts per million (CPM)") +
  scale_y_continuous(labels = scales::comma) +
  theme(
    legend.title = element_text(size = 16, face="bold"), 
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    strip.background = element_blank(), 
    strip.text.x = element_text(size = 14, face = "bold"))
DNA_barplot_virus_norm


##Create RNA phyloseq object
otu_table_RNA_virus <- RNA_virus_table %>% setNames(str_replace(colnames(.), "^RNA_(\\d+).*", "RNA\\1")) %>% otu_table(taxa_are_rows = T)
sample_data_RNA_virus <- RNA_metadata %>% select(-Name) %>% column_to_rownames("Sample") %>% sample_data()
tax_table_RNA_virus <- taxo_virus_all %>% as.matrix() %>% tax_table()

RNA_phylo_virus <- merge_phyloseq(otu_table_RNA_virus, sample_data_RNA_virus, tax_table_RNA_virus) %>% tax_fix()

RNA_phylo_virus_rel  <- transform_sample_counts(RNA_phylo_virus, function(x) x / sum(x))
#333 RNA vMAGs


RNA_colors <- c(
  Chrymotiviricetes = "#FFED6F",
  Alsuviricetes = "#CCEBC5", 
  Pisoniviricetes =  "#BC80BD",
  Caudoviricetes = "#8DD3C7",
  "Riboviria Realm" = "#FCCDE5",
  Arfiviricetes = "#FFFFB3",
  Resentoviricetes = "#B3DE69",
  Magsaviricetes = "#FDB462",
  Other = "lightgray")

#Plotting
RNA_barplot_virus <- RNA_phylo_virus_rel %>%
  ps_arrange(Succession) %>%
  merge_samples("Succession", fun=mean) %>% 
  comp_barplot(
    tax_level = "Class", n_taxa = 8,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9,
    palette = RNA_colors) + 
  theme_classic() +
  xlab("") +
  ylab("Relative abundance\n(RNA vMAGs)") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14, face = "bold"))
RNA_barplot_virus


##Normalize per sequencing depth
RNA_otu <- as(otu_table(RNA_phylo_virus), "matrix")
if (taxa_are_rows(RNA_phylo_virus)) {
  RNA_otu <- t(RNA_otu)
}

RNA_depths <- sample_data(RNA_phylo_virus)$Depth
names(RNA_depths) <- rownames(sample_data(RNA_phylo_virus))
RNA_depths <- RNA_depths[rownames(RNA_otu)]
RNA_otu_norm <- sweep(RNA_otu, 1, RNA_depths, "/") * 1e6

RNA_phylo_virus_norm <- phyloseq(otu_table(RNA_otu_norm, taxa_are_rows = FALSE), tax_table(RNA_phylo_virus), sample_data(RNA_phylo_virus))

RNA_norm_class <- tax_glom(RNA_phylo_virus_norm, taxrank = "Class")
RNA_plot_df <- psmelt(RNA_norm_class)

RNA_top_classes <- RNA_plot_df %>%
  group_by(Class) %>%
  summarise(mean_ab = mean(Abundance), .groups = "drop") %>%
  slice_max(mean_ab, n = 8) %>%
  pull(Class)

RNA_plot_df <- RNA_plot_df %>%
  mutate(Class = as.character(Class), 
         Class = ifelse(Class %in% RNA_top_classes, Class, "Other"))

RNA_plot_df_mean <- RNA_plot_df %>%
  group_by(Succession, Class) %>%
  summarise(Mean = mean(Abundance, na.rm = TRUE), .groups = "drop")

RNA_plot_df_mean$Class <- factor(RNA_plot_df_mean$Class, levels = c("Chrymotiviricetes", "Alsuviricetes", "Pisoniviricetes", "Caudoviricetes", "Riboviria Realm", "Arfiviricetes", "Resentoviricetes", "Magsaviricetes", "Other"))

RNA_plot_df_mean$Succession <- factor(RNA_plot_df_mean$Succession, levels = c("5","10","14","19","24","26","28","33",
                                                                    "38","40","42","47","54","61","67","75",
                                                                    "82","87","89","91","96","103"))

RNA_barplot_virus_norm <- ggplot(RNA_plot_df_mean, aes(x = Succession, y = Mean, fill = Class)) +
  geom_bar(stat = "identity", colour = "grey5", size = 0.1, width = 0.9, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = RNA_colors) +
  theme_classic() +
  xlab("") +
  ylab("Transcripts per million (TPM)") +
  xlab("Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14, face = "bold"))
RNA_barplot_virus_norm

FigureS5 <- ggarrange(DNA_barplot_virus, DNA_barplot_virus_norm, RNA_barplot_virus, RNA_barplot_virus_norm, ncol = 1, nrow = 4, labels = c("A", "B", "C", "D"), font.label = list(size = 18), align = "v", common.legend =  F, legend = "right")
FigureS5
ggsave("Figures/FigureS5.pdf", FigureS5, width = 18.2, height = 18.2)


###Relative abundance Phylum and Class

#DNA
Top_phyla_DNA_virus <- tax_glom(DNA_phylo_virus_rel, taxrank = "Phylum") %>% 
  ps_melt() %>% 
  group_by(Phylum) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
Top_class_DNA_virus <- tax_glom(DNA_phylo_virus_rel, taxrank = "Class") %>% 
  ps_melt() %>% 
  group_by(Class) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
DNA_prc <- DNA_norm_class %>%
  psmelt() %>%
  group_by(Sample) %>%
  summarise(total_CPM = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  summarise(min_prc = min(total_CPM)/10000, max_prc = max(total_CPM)/10000)

#RNA
Top_phyla_RNA_virus <- tax_glom(RNA_phylo_virus_rel, taxrank = "Phylum") %>% 
  ps_melt() %>% 
  group_by(Phylum) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
Top_class_RNA_virus <- tax_glom(RNA_phylo_virus_rel, taxrank = "Class") %>% 
  ps_melt() %>% 
  group_by(Class) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
RNA_prc <- RNA_norm_class %>%
  psmelt() %>%
  group_by(Sample) %>%
  summarise(total_CPM = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  summarise(min_prc = min(total_CPM)/10000, max_prc = max(total_CPM)/10000)
