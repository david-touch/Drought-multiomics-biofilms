

###Tree bacterial MAGs and quality assessment

#load libraries
library(tidyverse)
library(reshape2)
library(dplyr)
library(ggtree)
library(scales)
library(ggtreeExtra)
library(microViz)
library(ggnewscale)
library(ggpubr)
library(cowplot)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)


#Load data
bacteria_taxonomy <- read.delim("Data/gtdbtk.bac120.summary.tsv", sep="\t")
archaea_taxonomy <- read.delim("Data/gtdbtk.ar53.summary.tsv", sep="\t")
CheckM2_MAGs <- read.delim("Data/CheckM2.tsv", sep="\t")
MAGs_Galah <- read_csv("Data/MAGs_Galah.csv")
Euk_MAGs <- read_csv("Data/EukMAGs_Whokaryote_80.csv")
MAGs_euk_HQ <- read_csv("Data/MAGs_euk_HQ.csv")
Bac_tree	<- read.tree("Data/Bac_tree")
Rel_abundance_DNA <- read.delim("Data/DNA_coverage.tsv", sep="\t") %>% dplyr::rename(DNA_trimmean = DNA_other_cat_fwd.fq.gz.Trimmed.Mean)
Rel_abundance_RNA <- read.delim("Data/RNA_coverage_cat.tsv", sep="\t") %>% dplyr::rename(RNA_trimmean = RNA_other_cat_fwd.fq.gz.Trimmed.Mean)


#Merge files
Taxonomy <- rbind(bacteria_taxonomy, archaea_taxonomy) %>% dplyr::rename(Name = user_genome)
NumberBins <- n_distinct(Taxonomy$Name)
NumberBins
MAGs_info <- left_join(CheckM2_MAGs, Taxonomy)

#Filter interesting info
MAGs_info_HQ_bac <- MAGs_info %>% 
  filter(! (Name %in% Euk_MAGs$MAGs)) %>% 
  filter(Name %in% MAGs_Galah$MAGs) %>% 
  mutate(Quality = case_when(Completeness >= 90 & Contamination <= 5 ~ "High",
                             Completeness >= 70 & Contamination <= 10 ~ "Medium",
                             TRUE ~ "Low"),
         GC_Content = GC_Content * 100)
MAGs_info_HQ_bac_df <- MAGs_info_HQ_bac %>% filter(Quality == "High"| Quality == "Medium") %>% select(Name, Completeness, Contamination, Total_Contigs, Contig_N50, Genome_Size, GC_Content, Quality, classification)
write.csv(MAGs_info_HQ_bac_df, "Results/Supplementary_file_1.csv")

MAGs_all_bac <- MAGs_info_HQ_bac %>% select(Name, Completeness, Contamination, Contig_N50, Genome_Size, GC_Content, classification, Quality)
MAGs_HQ_bac <- MAGs_all_bac %>% filter(Quality == "High"| Quality == "Medium")

MAGs_HQ_bac_df <- MAGs_HQ_bac %>%  
  select(-Contig_N50) %>%
  mutate(Genome_Size = Genome_Size/1000000)
MAGs_HQ_bac_df$Quality <- factor(MAGs_HQ_bac_df$Quality, levels = c("High", "Medium"))
write.csv(MAGs_HQ_bac_df, "Results/MAGs_HQ_bac_df.csv")


##Get all HQ MAGs
MAGs_HQ_all <- MAGs_info %>% 
  filter(Name %in% MAGs_euk_HQ$MAGs | Name %in% MAGs_HQ_bac$Name) %>%
  select(Name) %>% dplyr::rename(MAGs = Name)
write.csv(MAGs_HQ_all, "Results/MAGs_HQ_all.csv", row.names = FALSE)

p1 <- MAGs_HQ_bac_df %>% 
  ggboxplot(
    x = "Quality",
    y = c("Completeness", "GC_Content", "Contamination", "Genome_Size"),
    fill = "Quality",
    palette =c("forestgreen", "greenyellow"),
    title = "",
    xlab = "",
    ylab = c(
      "Completeness (%)",
      "GC content (%)",
      "Contamination (%)",
      "Genome size (Mbp)"),
    ggtheme = theme(panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size=14, face = "bold"),
                    axis.line = element_line(color = "black"),
                    axis.title = element_text(size = 14, color = "black"),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 0),
                    axis.ticks.y = element_blank(),
                    plot.background = element_blank()))

p1$Completeness <- p1$Completeness + rotate()
p1$GC_Content <- p1$GC_Content + rotate()
p1$Contamination <- p1$Contamination + rotate() + scale_y_continuous(limits = c(0,10))
p1$Genome_Size <- p1$Genome_Size + rotate() + scale_y_continuous(limits = c(0,10))

FigureS1 <- ggarrange(plotlist = p1, common.legend = TRUE)
FigureS1
ggsave("Figures/FigureS1.pdf", FigureS1, width = 9.1, height = 9.1)


#Info
Size_mean <- mean(MAGs_HQ_bac_df$Genome_Size)
Size_mean
Size_sd <- sd(MAGs_HQ_bac_df$Genome_Size)
Size_sd

GC_mean <- mean(MAGs_HQ_bac_df$GC_Content)
GC_mean
GC_sd <- sd(MAGs_HQ_bac_df$GC_Content)
GC_sd

HQ_taxonomy <- Taxonomy %>% select(Name, classification) %>% 
  filter(Name %in% MAGs_HQ_bac$Name) %>% 
  dplyr::rename(MAG = Name) %>% 
  separate(classification, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
HQ_taxonomy_long <- HQ_taxonomy %>% 
  mutate(Phylum = case_when(Phylum == "p__Pseudomonadota" ~ Class, TRUE ~ Phylum)) %>% 
  pivot_longer(cols = !MAG, names_to = "Taxa", values_to = "Name") %>% 
  mutate(Name = gsub(".__", "", Name)) %>% 
  mutate(Name = gsub("_.", "", Name)) %>% 
  mutate(Status = if_else(Name == "", "UnKnown", "Known"))
HQ_taxonomy_count <- HQ_taxonomy_long %>%
  select(-MAG) %>%
  group_by(Taxa, Status) %>%
  summarise(
    sum = n(),
    perc = sum / 117,
    label = percent(perc, accuracy = 0.1))

#91.5% of the HQ MAGs are identified at the genus level, 100% are at the family level.


HQ_phyla	<- HQ_taxonomy_long %>%
  filter(Taxa == "Phylum") %>% 
  mutate(MAGS = MAG)

allPhyla <- HQ_phyla %>% 
  group_by(Name) %>%
  summarise(total	= n()) %>%
  arrange(desc(total))
allPhyla

list_phyla <- allPhyla$Name
a_phyla	<- split(HQ_phyla$MAG, HQ_phyla$Name)

mostAbundantFamily	<- HQ_taxonomy_long %>%
  filter(Taxa == "Family") %>% 
  group_by(Name) %>%
  summarise(total	= n()) %>%
  arrange(desc(total)) %>%
  dplyr::slice(1:9)
mostAbundantFamily

list_family <- mostAbundantFamily$Name

HQ_family	<- HQ_taxonomy_long %>%
  filter(Taxa == "Family") %>% 
  select(-Status) %>% 
  mutate(Name	= if_else(Name %in% list_family, Name, "Other"))

a_family	<- split(HQ_family$MAG, HQ_family$Name)
a_family
list_family[10] <- "Other"
a_family_ordered <- a_family[match(list_family, names(a_family))]

tree_family	<- groupOTU(Bac_tree, a_family_ordered)

family_color_map <- distinct_palette(n = 9, pal = "kelly")
family_color_map <- setNames(family_color_map[1:length(list_family)], list_family)
family_color_map <- c("lightgray", family_color_map)
family_color_map["Burkholderiaceae"] <- "#875692"
family_color_map["Methylophilaceae"] <- "#f3c300"

names(family_color_map) <- c("A", names(family_color_map)[-1])
ColorFam <- family_color_map[order(names(family_color_map))]
ColorFam <- unname(ColorFam)

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

Rel_abundance_percent <- Rel_abundance_DNA %>% column_to_rownames("Genome") %>% 
  mutate(Rel_abundance_DNA = (DNA_trimmean/sum(DNA_trimmean))*100) %>%
  rownames_to_column("Genome") %>% 
  left_join(Rel_abundance_RNA) %>% 
  mutate(Rel_abundance_RNA = (RNA_trimmean/sum(RNA_trimmean))*100) %>% 
  mutate(FC = (Rel_abundance_RNA - Rel_abundance_DNA)/Rel_abundance_DNA) %>% 
  select(-DNA_trimmean, -RNA_trimmean)

HQ_phyla <- HQ_phyla %>% left_join(Rel_abundance_percent, by = join_by(MAG == Genome))
HQ_family <- HQ_family %>% left_join(Rel_abundance_percent, by = join_by(MAG == Genome))

Rel_abundance_percent_filt <- Rel_abundance_percent %>% filter(Genome %in% HQ_phyla$MAG)

DNA_ring <- Rel_abundance_percent_filt %>% column_to_rownames("Genome") %>% select(Rel_abundance_DNA)
RNA_ring <- Rel_abundance_percent_filt %>% column_to_rownames("Genome") %>% select(Rel_abundance_RNA)
FC_ring <- Rel_abundance_percent_filt %>% column_to_rownames("Genome") %>% select(FC)

Rel_ab_Gammaproteobacteria <- HQ_phyla %>% filter(Name == "Gammaproteobacteria") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Gammaproteobacteria
Rel_ab_Alphaproteobacteria <- HQ_phyla %>% filter(Name == "Alphaproteobacteria") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Alphaproteobacteria
Rel_ab_Bacteroidota <- HQ_phyla %>% filter(Name == "Bacteroidota") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Bacteroidota

Rel_ab_Burkholderiaceae <- HQ_family %>% filter(Name == "Burkholderiaceae") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Burkholderiaceae
Rel_ab_Methylophilaceae <- HQ_family %>% filter(Name == "Methylophilaceae") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Methylophilaceae
Rel_ab_Sphingomonadaceae <- HQ_family %>% filter(Name == "Sphingomonadaceae") %>% summarise(sum = sum(Rel_abundance_DNA, na.rm = TRUE)) %>% pull(sum)
Rel_ab_Sphingomonadaceae

f1	<-
  ggtree(tree_family, layout = 'circular', branch.length="none", aes(color = group)) +
  geom_tree() +
  theme_tree() + 
  theme(legend.position = "right") +
  scale_color_manual(values	= ColorFam,
                     na.value = "transparent",
                     guide = "none")
f1

f2 <- f1 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = HQ_family,
    pwidth	= 0.03,
    geom = geom_tile,
    mapping = aes(y = MAG, fill = Name, x = 2), stat = "identity") +
  scale_fill_manual(values	= family_color_map) +
  scale_y_continuous(expand = c(-0.0025,0)) +
  guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))
f2

f3 <- f2 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = HQ_phyla,
    pwidth	= 0.03,
    geom = geom_tile,
    mapping = aes(y = MAG, fill = factor(Name, levels = names(ColorPhyla)), x = 1),
    offset = 0.3,
    stat = "identity") +
  scale_fill_manual(values	= ColorPhyla) +
  scale_y_continuous(expand = c(-0.0025,0)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=14, face = "bold"),
        plot.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(title = "Phylum", reverse = FALSE))
f3
leg3 <- get_legend(f3)
leg_3 <- as_ggplot(leg3)

f3 <- f3 +
  guides(fill = "none")

f4 <- f3 +
  new_scale_fill() +
  geom_fruit(
    data = HQ_phyla,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = MAG, fill = Rel_abundance_DNA, x = 1),
    offset = 0.01,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none") +
  new_scale_colour() +
  new_scale_fill()
f4  

f5 <- 
  gheatmap(
    f4,
    DNA_ring,
    offset = 12,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "white", mid = "#a6cee3", high = "#1F78B4", midpoint = 2, limits = c(0,6), breaks = c(0, 2, 4, 6)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_colorbar(title = "gDNA coverage (%)", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f5  
leg5 <- get_legend(f5)
leg_5 <- as_ggplot(leg5) + theme(plot.margin = margin(0, 0, 0, 0))

f5 <- f5 + guides(fill = "none")

f6 <- f5 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = HQ_phyla,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = MAG, fill = Rel_abundance_RNA, x = 1),
    offset = 0,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none") +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  new_scale_colour() +
  new_scale_fill()
f6

f7 <- 
  gheatmap(
    f6,
    RNA_ring,
    offset = 14,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "white", mid = "#fb9a99", high = "#E31A1C", midpoint = 3.5, limits = c(0, 8), breaks = c(0, 2, 4, 6, 8))+
  theme(plot.margin = margin(0, 0, 0, 0)) +
  labs(fill = "Relative\nabundance (%)") +
  guides(fill = guide_colorbar(title = "mRNA coverage (%)", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f7
leg7 <- get_legend(f7)
leg_7 <- as_ggplot(leg7) + theme(plot.margin = margin(0, 0, 0, 0))

f7 <- f7 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

###
f8 <- f7 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = HQ_phyla,
    pwidth	= 0.03,
    geom = geom_bar,
    mapping = aes(y = MAG, fill = FC, x = 1),
    offset = 0,
    stat = "identity", 
    alpha = 0) +
  guides(fill = "none")
f8

f9 <- 
  gheatmap(
    f8,
    FC_ring,
    offset = 16,
    width = 0.05,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_gradient2(low = "purple", mid = "white", high = "darkgreen", midpoint = 0, na.value = "white", limits = c(-1, 6), breaks = c(-1, 2, 4, 6)) +
  guides(fill = guide_colorbar(title = "Fold change", 
                               barwidth = 6, 
                               barheight = 1, 
                               direction = "horizontal",
                               title.position = "top"))
f9
leg9 <- get_legend(f9)
leg_9 <- as_ggplot(leg9) + theme(plot.margin = margin(0, 0, 0, 0))

f9 <- f9 + guides(fill = "none") + theme(plot.margin = margin(0, 0, 0, 0))

legends <- ggarrange(leg_3, leg_5, leg_7, leg_9, ncol = 1, heights = c(0.6, 0.2, 0.2, 0.2))
legends

f10 <- ggarrange(f9, legends, heights = c(1, 0.5))
f10

ggsave("Figures/Figure1.pdf", f10)
