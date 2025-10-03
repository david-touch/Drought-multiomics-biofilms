

###Tree eukaryotic MAGs

#load libraries
library(tidyverse)
library(ggtree)
library(scales)
library(ggtreeExtra)
library(ggnewscale)
library(ggpubr)
library(ape)
library(microViz)
library(cowplot)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data
Euk_tree	<- read.tree("Data/Euk_tree")
Euk_MAGs <- read_csv("Data/EukMAGs_Whokaryote_80.csv")
Rel_abundance_DNA <- read.delim("Data/DNA_coverage.tsv", sep="\t") %>% dplyr::rename(DNA_trimmean = DNA_other_cat_fwd.fq.gz.Trimmed.Mean)
Rel_abundance_RNA <- read.delim("Data/RNA_coverage_cat.tsv", sep="\t") %>% dplyr::rename(RNA_trimmean = RNA_other_cat_fwd.fq.gz.Trimmed.Mean)
Euk_ref <- read.delim("Data/euk_metadata_REF.tsv", sep="\t")
Euk_input <- read.delim("Data/euk_input_metadata.tsv", sep="\t")


#Modify dataframe
Euk_input <- Euk_input %>% mutate(File.Name = gsub("\\_results.fas$", "", File.Name),
                                  File.Name = gsub("\\.fas$", "", File.Name))
Euk_input_mod <- Euk_input %>% select(File.Name, Unique.ID, Higher.Taxonomy, Lower.Taxonomy) %>% 
  mutate(Lower.Taxonomy = gsub("Ochrophyta-Bolidophyceae", "Ochrophyta-Bacillariophyceae", Lower.Taxonomy),
         Type = "MAG")
Euk_df <- Euk_ref %>% select(Unique.ID, Higher.Taxonomy, Lower.Taxonomy) %>% mutate(Type = "ref", File.Name = "")
Euk_df <- rbind(Euk_df, Euk_input_mod)

#Filter HQ euk
MAGs_euk_HQ <- Euk_input %>% select(File.Name, Higher.Taxonomy, Lower.Taxonomy) %>% dplyr::rename(MAGs = File.Name, Phylum = Higher.Taxonomy, Family = Lower.Taxonomy) %>% 
  filter(MAGs != "Diatom1" & MAGs != "Diatom2" & MAGs != "Hydrurus")
write.csv(MAGs_euk_HQ, "Results/MAGs_euk_HQ.csv", row.names = FALSE)

#Filter 3 added reference genomes
Euk_tree_filt <- drop.tip(Euk_tree, c("R1", "R2", "R3"))
Euk_df_filt <- Euk_df %>% filter(!Unique.ID %in% c("R1", "R2", "R3"))

#Plot tree
taxToPlot <- Euk_df_filt %>%
  filter(Type == "MAG") %>%
  filter(Unique.ID %in% Euk_tree_filt$tip.label)

Euk_df_filt <- Euk_df_filt %>%
  mutate(taxa_good = if_else(Lower.Taxonomy %in% taxToPlot$Lower.Taxonomy, Lower.Taxonomy, "Other"))

tax_temp <- Euk_df_filt %>%
  filter(Unique.ID %in% Euk_tree$tip.label)

tree_meta	<- split(Euk_df_filt$Unique.ID, Euk_df_filt$taxa_good)
Euk_tree_lower	<- groupOTU(Euk_tree_filt, tree_meta)

EukColorClass <- c("#b3446c", "#8db600", "lightgrey")

mostabundantHigher <- Euk_df %>% 
  group_by(Higher.Taxonomy) %>%
  summarise(total	= n()) %>%
  arrange(desc(total)) %>%
  dplyr::slice(1:10)

list_higher <- mostabundantHigher$Higher.Taxonomy
list_higher[11] <- "Other"

df_higher	<- Euk_df %>%
  select(-Lower.Taxonomy, -Type) %>% 
  mutate(Higher.Taxonomy	= if_else(Higher.Taxonomy %in% list_higher, Higher.Taxonomy, "Other")) %>% 
  dplyr::rename(HigherTaxonomyTemp = Higher.Taxonomy)
df_higher <- as_tibble(df_higher)

EukColorHigher <- distinct_palette(n = 10, pal = "greenArmytage")
EukColorHigher <- setNames(EukColorHigher[1:length(list_higher)], list_higher)
EukColorHigher["Stramenopiles"] <- "#33A02C"
EukColorHigher["Rhizaria"] <- "#F0A3FF"
EukColorHigher["Rhodophyta"] <- "#C20088"
EukColorHigher["Obazoa"] <- "#1ff8ff"

Euk_rel_abundance_percent <- Rel_abundance_DNA %>% column_to_rownames("Genome") %>% 
  mutate(Rel_abundance_DNA = (DNA_trimmean/sum(DNA_trimmean))*100) %>%
  rownames_to_column("Genome") %>% 
  left_join(Rel_abundance_RNA) %>% 
  mutate(Rel_abundance_RNA = (RNA_trimmean/sum(RNA_trimmean))*100) %>% 
  mutate(FC = (Rel_abundance_RNA - Rel_abundance_DNA)/Rel_abundance_DNA) %>% 
  select(-DNA_trimmean, -RNA_trimmean) %>% 
  filter(Genome %in% MAGs_euk_HQ$MAGs)

df_higher <- df_higher %>% left_join(Euk_rel_abundance_percent, by = join_by(File.Name == Genome))
Euk_rel_abundance_percent_DNA <- df_higher %>% column_to_rownames("Unique.ID") %>% select(Rel_abundance_DNA)
Euk_rel_abundance_percent_RNA <- df_higher %>% column_to_rownames("Unique.ID") %>% select(Rel_abundance_RNA)
Euk_rel_abundance_percent_FC <- df_higher %>% column_to_rownames("Unique.ID") %>% select(FC)

Rel_ab_Stramenopiles <- df_higher %>% filter(HigherTaxonomyTemp == "Stramenopiles") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Stramenopiles

f1	<-
  ggtree(Euk_tree_lower, layout = 'circular', branch.length="none", aes(color = group)) +
  geom_tree() +
  theme_tree() + 
  theme(legend.position = "right",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_color_manual(values	= EukColorClass,
                     na.value = "transparent",
                     guide = "none")
f1

f2 <- f1 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = Euk_df_filt,
    pwidth	= 0.03,
    geom = geom_tile,
    mapping = aes(y = Unique.ID, fill = taxa_good, x = 2), stat = "identity") +
  scale_fill_manual(values	= EukColorClass) +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  scale_y_continuous(expand = c(-0.0025,0)) +
  guides(fill = "none")
f2

f3 <- f2 %<+% tax_temp +
  geom_tippoint(aes(color = Type, size = Type, shape = Type)) +
  scale_color_manual(values = c("#004488", "white")) +
  scale_size_manual(values = c(3, 0)) +
  scale_shape_manual(values = c(1, 16)) +
  guides(color = "none", size = "none", shape = "none") +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  new_scale_colour() +
  new_scale_fill()
f3

f4 <- f3 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = df_higher,
    pwidth	= 0.03,
    geom = geom_tile,
    mapping = aes(y = Unique.ID, fill = factor(HigherTaxonomyTemp, levels = names(EukColorHigher)), x = 1),
    offset = 0.3,
    stat = "identity") +
  scale_fill_manual(values	= EukColorHigher) +
  scale_y_continuous(expand = c(-0.0025,0)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=14, face = "bold"),
        plot.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(title = "Phylum", reverse = FALSE))
f4

leg4 <- get_legend(f4)
leg_4 <- as_ggplot(leg4)

f4 <- f4 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

f5 <- f4 +
  new_scale_colour() +
  new_scale_fill() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  geom_fruit(
    data = df_higher,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = Unique.ID, fill = Rel_abundance_DNA, x = 1),
    offset = 0,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none") +
  new_scale_colour() +
  new_scale_fill()
f5

f6 <- 
  gheatmap(
    f5,
    Euk_rel_abundance_percent_DNA,
    offset = 10,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "white", mid = "#a6cee3", high = "#1F78B4", midpoint = 2, na.value = "transparent", limits = c(0,6), breaks = c(0, 2, 4, 6)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_colorbar(title = "gDNA coverage (%)", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f6  

leg6 <- get_legend(f6)
leg_6 <- as_ggplot(leg6)

f6 <- f6 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

f7 <- f6 +
  new_scale_colour() +
  new_scale_fill() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  geom_fruit(
    data = df_higher,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = Unique.ID, fill = Rel_abundance_RNA, x = 1),
    offset = 0,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none") +
  new_scale_colour() +
  new_scale_fill()
f7

f8 <- 
  gheatmap(
    f7,
    Euk_rel_abundance_percent_RNA,
    offset = 12,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "white", mid = "#fb9a99", high = "#E31A1C", midpoint = 3.5, na.value = "transparent", limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  labs(fill = "Relative\nabundance (%)") +
  guides(fill = guide_colorbar(title = "mRNA coverage (%)", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f8
leg8 <- get_legend(f8)
leg_8 <- as_ggplot(leg8) + theme(plot.margin = margin(0, 0, 0, 0))

f8 <- f8 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

legends_euk <- ggarrange(leg_4, leg_6, leg_8, ncol = 1, heights = c(0.6, 0.2, 0.2))
legends_euk

f9 <- f8 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = df_higher,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = Unique.ID, fill = FC, x = 1),
    offset = 0,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none")
f9

f10 <- 
  gheatmap(
    f9,
    Euk_rel_abundance_percent_FC,
    offset = 14,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "purple", mid = "white", high = "darkgreen", midpoint = 0, na.value = "transparent",  limits = c(-1, 6), breaks = c(-1, 2, 4, 6)) +
  guides(fill = guide_colorbar(title = "Fold change", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f10
leg10 <- get_legend(f10)
leg_10 <- as_ggplot(leg10) + theme(plot.margin = margin(0, 0, 0, 0))

f11 <- f10 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

legends <- ggarrange(leg_4, leg_6, leg_8, leg_10, ncol = 1, heights = c(0.6, 0.2, 0.2, 0.2))
legends

f12 <- ggarrange(f11, legends, heights = c(1, 0.5))
f12

ggsave("Figures/FigureS2.pdf", f12)
