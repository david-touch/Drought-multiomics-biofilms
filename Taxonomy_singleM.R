

###Taxonomic changes on whole metagenome and whole metatranscpritomic 

#Load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)
library(purrr)
library(broom)
library(grid)
library(cowplot)
library(ggh4x)


#Choose directory
setwd("~/Drought-multiomics-biofilms")
set.seed(19950930)

#Load and prepare data
DNA_metadata <- read_csv("Data/DNA_metadata.csv") %>% data.frame()
singlem_DNA <- read.delim("Data/DNA_otu_table_singlem.tsv", sep="\t")

singlem_DNA_combined <- singlem_DNA %>% 
  group_by(sample, taxonomy) %>% 
  summarise(SumCov = sum(coverage)) %>% 
  filter(taxonomy != "Root") %>% 
  pivot_wider(names_from = sample, values_from = SumCov)

singlem_DNA_combined_mat <- singlem_DNA_combined %>% 
  column_to_rownames("taxonomy") %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  as.matrix()

#Create a phyloseq object
Phototroph_euk <- c("Cryptophyceae", "Glaucocystophyceae", "Haptista", "Rhodophyta", "Sar", "Viridiplantae")

otu_table_DNA_singlem <- singlem_DNA_combined_mat %>% otu_table(taxa_are_rows = T)
sample_data_DNA_singlem <- DNA_metadata %>% mutate(Name = str_remove(Name, "\\.fq\\.gz\\.Trimmed\\.Mean$")) %>% column_to_rownames("Name") %>% mutate(ID = Sample) %>% sample_data()
tax_table_DNA_singlem <- singlem_DNA_combined %>% select(taxonomy) %>% mutate(taxo = taxonomy) %>% column_to_rownames("taxonomy") %>% 
  separate(taxo, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = "; ") %>% mutate(across(Domain:Genus, ~str_remove(., "^[a-z]__"))) %>%
  select(-Root) %>% mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Eukaryota" & Phylum %in% Phototroph_euk ~ "Phototrophic eukaryotes",
                                              Domain == "Eukaryota" & !Phylum %in% Phototroph_euk ~ "Non-phototrophic eukaryotes",
                                              TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Bacteria" & is.na(Phylum) ~ "Unclassified", TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Archaea" ~ "Archaea Domain", TRUE ~ Phylum),
                           Phylum = case_when(Phylum == "Patescibacteriota" ~ "Patescibacteria", TRUE ~ Phylum)) %>% 
  as.matrix() %>% tax_table()

DNA_phylo_singlem <- merge_phyloseq(otu_table_DNA_singlem, sample_data_DNA_singlem, tax_table_DNA_singlem) %>% tax_fix()

DNA_phylo_singlem_rel  <- transform_sample_counts(DNA_phylo_singlem, function(x) x / sum(x))

ColorPhyla <- c(
  Gammaproteobacteria = "#A6CEE3",
  Alphaproteobacteria = "#1F78B4",
  Bacteroidota = "#FFFF99",
  Cyanobacteriota = "#B2DF8A", 
  "Phototrophic eukaryotes" = "#33A02C",
  "Non-phototrophic eukaryotes" = "#1ff8ff",
  Patescibacteria = "#FB9A99",
  Actinomycetota = "#E31A1C",
  Verrucomicrobiota = "#FDBF6F", 
  Planctomycetota = "#FF7F00",  
  Chloroflexota = "#CAB2D6",
  Spirochaetota = "#6A3D9A", 
  Bdellovibrionota = "#1B9E77",
  Fibrobacterota = "#B15928", 
  Myxococcota = "#7570B3",
  Deinococcota = "#D95f02",
  Bacillota = "#E7298A",
  Acidobacteriota = "#66A61E",
  "Archaea Domain" = "#E6AB02",
  Other = "lightgray", 
  Unclassified = "#5C5C5C")

Phyla_order <- c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidota", "Cyanobacteriota", "Phototrophic eukaryotes", "Non-phototrophic eukaryotes" , "Patescibacteria", "Actinomycetota", "Verrucomicrobiota", "Planctomycetota", "Chloroflexota", "Spirochaetota", "Bdellovibrionota", "Fibrobacterota", "Obazoa", "Myxococcota", "Deinococcota", "Bacillota", "Acidobacteriota", "Archaea", "Other")

Label_droughts = as_labeller(c("1" = "Drought 1 (6h)", "2" = "Drought 2 (24h)", "3" = "Drought 3 (24h)"))

Label_samples_DNA <- c("S4" = "Before", "S5" = "After",
                       "S8" = "Before", "S9" = "After",
                       "S17" = "Before", "S18" = "After")

#Plotting
DNA_barplot_singlem_all <- DNA_phylo_singlem_rel %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = TRUE,
    sample_order = "default",
    bar_width = 0.9,
    palette = ColorPhyla) + 
  theme_classic() +
  xlab("Days of growth")+
  ylab("Relative abundance") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"))
DNA_barplot_singlem_all

ggsave("Figures/FigureS3.pdf", DNA_barplot_singlem_all, width = 18.2, height = 9.1)


#Relative abundance Kingdom and Phyla
Kingdom_singlem <- tax_glom(DNA_phylo_singlem_rel, taxrank = "Domain") %>% 
  ps_melt() %>% 
  group_by(Domain) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))
Phyla_singlem <- tax_glom(DNA_phylo_singlem_rel, taxrank = "Phylum") %>% 
  ps_melt() %>% 
  group_by(Phylum) %>% 
  summarise(av_abund = mean(Abundance), sd_abund = sd(Abundance))%>% 
  arrange(desc(av_abund))


##Whole metatranscriptome taxonomy
RNA_metadata <- read.csv("Data/RNA_metadata.csv") %>% filter(Sample != "RNANEG") %>% 
  filter(Sequence == "sorted" & Target == "MAGs") %>% select(-Sequence, -Target)

RNA_metadata$Sample <- factor(RNA_metadata$Sample)
RNA_metadata$Time <- factor(RNA_metadata$Time)
RNA_metadata$Drought <- factor(RNA_metadata$Drought, levels = c("B1","A1", "B2", "A2", "B3", "A3"))
RNA_metadata$State <- factor(RNA_metadata$State, levels = c("Before","After"))
RNA_metadata$Period <- as.factor(RNA_metadata$Period)

singlem_RNA <- read.delim("Data/RNA_otu_table_singlem.tsv", sep="\t") %>% 
  filter(sample != "RNA_other_cat_fwd")

singlem_RNA_combined <- singlem_RNA %>% 
  group_by(sample, taxonomy) %>% 
  summarise(SumCov = sum(coverage)) %>% 
  filter(taxonomy != "Root") %>% 
  pivot_wider(names_from = sample, values_from = SumCov)

singlem_RNA_combined_mat <- singlem_RNA_combined %>% 
  column_to_rownames("taxonomy") %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  rename_with(~ str_remove(., "_other_paired_fwd")) %>% 
  as.matrix()

#Create a phyloseq object
otu_table_RNA_singlem <- singlem_RNA_combined_mat %>% otu_table(taxa_are_rows = T)
sample_data_RNA_singlem <- RNA_metadata %>% mutate(Name = str_remove(Name, "_sorted_MAGs$")) %>% column_to_rownames("Name") %>% mutate(ID = Sample) %>% sample_data()
tax_table_RNA_singlem <- singlem_RNA_combined %>% select(taxonomy) %>% mutate(taxo = taxonomy) %>% column_to_rownames("taxonomy") %>% 
  separate(taxo, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = "; ") %>% mutate(across(Domain:Genus, ~str_remove(., "^[a-z]__"))) %>%
  select(-Root) %>% mutate(Phylum = case_when(Phylum == "Pseudomonadota" ~ Class, TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Eukaryota" & Phylum %in% Phototroph_euk ~ "Phototrophic eukaryotes",
                                              Domain == "Eukaryota" & !Phylum %in% Phototroph_euk ~ "Non-phototrophic eukaryotes",
                                              TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Bacteria" & is.na(Phylum) ~ "Other", TRUE ~ Phylum),
                           Phylum = case_when(Domain == "Archaea" ~ "Archaea Domain", TRUE ~ Phylum),
                           Phylum = case_when(Phylum == "Patescibacteriota" ~ "Patescibacteria", TRUE ~ Phylum),
                           Phylum = case_when(Phylum == "Bacteroidota_A" ~ "Bacteroidota", TRUE ~ Phylum)) %>% 
  as.matrix() %>% tax_table()

RNA_phylo_singlem <- merge_phyloseq(otu_table_RNA_singlem, sample_data_RNA_singlem, tax_table_RNA_singlem) %>% tax_fix()

RNA_phylo_singlem_rel  <- transform_sample_counts(RNA_phylo_singlem, function(x) x / sum(x))

Label_samples_RNA <- c("B1" = "Before", "A1" = "After",
                       "B2" = "Before", "A2" = "After",
                       "B3" = "Before", "A3" = "After")

#Plotting
RNA_barplot_droughts_singlem <- RNA_phylo_singlem_rel %>%
  ps_filter(!is.na(Period)) %>% 
  #merge_samples("Drought", fun = "mean") %>% 
  ps_arrange(State) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 20,
    tax_order = Phyla_order,
    x = "SAMPLE",
    bar_outline_colour = "grey5",
    merge_other = TRUE,
    sample_order = "default",
    bar_width = 0.9,
    palette = ColorPhyla) + 
  theme_classic() +
  xlab("")+
  ylab("Relative abundance") +
  scale_x_discrete(labels = Label_samples_RNA) +
  facet_nested(~ Period + State, scales = "free_x", space = "free_x", nest_line = element_line(linetype = 1), labeller = labeller(Period = Label_droughts)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14),
        axis.text.y =  element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold"))
RNA_barplot_droughts_singlem

FigureS7 <- RNA_barplot_droughts_singlem
FigureS7
ggsave("Figures/FigureS7.pdf", FigureS7, width = 18.2, height = 9.1)


#Ratio phototrophs 
tax_table_DNA_singlem_photo <- singlem_DNA_combined %>% select(taxonomy) %>% mutate(taxo = taxonomy) %>% column_to_rownames("taxonomy") %>% 
  separate(taxo, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = "; ") %>% mutate(across(Domain:Genus, ~str_remove(., "^[a-z]__"))) %>%
  select(-Root) %>% filter(Phylum == "Cryptophyceae" | Phylum == "Glaucocystophyceae" | Phylum == "Haptista" | Phylum == "Rhodophyta" | Phylum == "Sar" | Phylum == "Viridiplantae" | Phylum == "Cyanobacteriota" | Phylum == "Chloroflexota") %>% 
  mutate(Phylum = case_when(Domain == "Eukaryota" ~ "Phototrophic eukaryotes", TRUE ~ Phylum)) %>% 
  as.matrix() %>% tax_table()

DNA_phylo_singlem_photo <- merge_phyloseq(otu_table_DNA_singlem, sample_data_DNA_singlem, tax_table_DNA_singlem_photo) %>% tax_fix()

Phototrophs <- DNA_phylo_singlem_photo %>% 
  tax_select(tax_list = c("Cyanobacteriota", "Phototrophic eukaryotes")) %>% 
  tax_glom("Phylum")

ChlA <- DNA_metadata %>% 
  mutate(ChlA = rowMeans(select(., ChlA1, ChlA2, ChlA3))) %>% 
  select(Sample, Succession, Time, ChlA)

Phototrophs_ratio <- data.frame(otu_table(Phototrophs)) %>% 
  rownames_to_column("Phylum") %>% 
  mutate(Phylum = case_when(Phylum == "Root; d__Bacteria; p__Cyanobacteriota; c__Cyanobacteriia; o__Cyanobacteriales; f__Microcoleaceae; g__Microcoleus" ~ "Cyanobacteriota", TRUE ~ "Eukaryota")) %>% 
  {setNames(as.numeric(.[.$Phylum == "Eukaryota", -1]) / as.numeric(.[.$Phylum == "Cyanobacteriota", -1]), names(.[-1]))} %>% 
  enframe(name = "Sample", value = "Ratio") %>% 
  mutate(Sample = paste0("S", str_extract(Sample, "(?<=Sentinel_)\\d+")))

Phototrophs_ratio <- Phototrophs_ratio %>% left_join(ChlA)

#Plot ratio
Phototrophs_ratio$Sample <- factor(Phototrophs_ratio$Sample, levels = Phototrophs_ratio$Sample[order(Phototrophs_ratio$Time)])

max_ratio <- max(Phototrophs_ratio$Ratio, na.rm = TRUE)
max_chla <- max(Phototrophs_ratio$ChlA, na.rm = TRUE)
scale_factor <- max_ratio / max_chla

Phototrophs_ratio$ChlA_scaled <- Phototrophs_ratio$ChlA * scale_factor

ratio_plot <- Phototrophs_ratio %>% 
  ggplot(aes(x = Succession)) +
  geom_point(aes(y = Ratio), color = "#6A3D9A", size = 3, alpha = 0.5) +
  geom_smooth(aes(y = Ratio), se = FALSE, color = "#6A3D9A") +
  geom_point(aes(y = ChlA_scaled), color = "#006400", shape = 17, size = 3, alpha = 0.5) +
  geom_smooth(aes(y = ChlA_scaled), se = FALSE, color = "#006400", linetype = "dashed", linewidth = 1.5) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 105)) +
  scale_y_log10(name = "Ratio Phototrophic eukaryotes:Cyanobacteria",
                sec.axis = sec_axis(~ . / scale_factor, name = bquote(bold('Chlorophyll-a concentration'~(µg/cm^2))))) +
  theme_classic() +
  xlab("Succession (days)") +
  theme(
    legend.title = element_text(size = 16, face="bold"), 
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold", margin = margin(r = 10)),
    axis.text.x =  element_text(size = 14),
    axis.text.y =  element_text(size = 12),
    strip.background = element_blank(), 
    strip.text.x = element_text(size = 14, face = "bold"))
ratio_plot

g <- ggplotGrob(ratio_plot)
g$grobs[[which(g$layout$name == "ylab-l")]]$children[[1]]$gp$col <- "#6A3D9A"
g$grobs[[which(g$layout$name == "ylab-r")]]$children[[1]]$gp$col <- "#006400"
grid.newpage()
grid.draw(g)                
Photo_ratio <- ggdraw(g)                
Photo_ratio

ggsave("Figures/Figure6.pdf", Photo_ratio, width = 9.1, height = 7)
