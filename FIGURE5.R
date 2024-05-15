############################################################################################################################
#### CODE USED TO GENERATE FIGURE 5
############################################################################################################################

# FIGURE 5A

library(pheatmap)
matrix <- read.delim("RPKMcounts.txt", header = TRUE, sep = "\t", row.names = 1)
matrix <- na.omit(matrix)
transformed_matrix <- log10(matrix+0.01)
my_sample_col <- read.table("metadata.txt", header = TRUE, sep = "\t",row.names = 1)

family_sums <- rowSums(transformed_matrix)
top10_families <- names(sort(family_sums, decreasing = TRUE))[1:10]
matrix_top10 <- transformed_matrix[top10_families, ]


ALR_heatmap <- pheatmap(
  matrix_top10,
  color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = my_sample_col,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  fontsize = 20,
  fontsize_row = 30,
  legend = TRUE,
  display_numbers = F,
  annotation_legend = TRUE, show_colnames = FALSE, show_rownames = T
)
pdf("ALR_heatmap_prevhap.pdf",width=20,height=14);
ALR_heatmap
dev.off()

# FIGURE 5B - Fisher HAP/no HAP contigs

#Use Fisher test to identify discriminant contigs in HAP and no HAP signature 6 days before the HAP onset.
#Input table : Presence absence counts

otu_data <- read.delim("RPKM_NEW_WTA_deconta_B_50_01.txt", header = TRUE)

for (n in 1:nrow(otu_data)){
  otu_data_2<-otu_data[n,-1]
  otu_data_3<-rbind(otu_data_2,apply(otu_data[-n,-1],2,FUN=sum))
  otu_data$p.value[n]<-fisher.test(otu_data_3)$p.value
}

library(ggplot2)

data <- read.delim("fisher_output.txt")

HAP_noHAP_signature_contigs <- ggplot(data, aes(reorder(OTU, PVAL), PVAL, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.8, size = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey")) +
  theme(strip.placement = "outside",
                     strip.text.y = element_text(angle = 0)) +
  labs(title = "", x = "Significant contigs", y = "Fisher test -log10(P-Value)")+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text = element_text(size=30),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        plot.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = "bottom") + facet_wrap(~GROUP)


pdf("HAP_noHAP_signature_Fisher_PREVHAP.pdf",width=16,height=11);
HAP_noHAP_signature_contigs
dev.off()


# FIGURE 5C - Lefse HAP contigs

#Use LEfSe to identify discriminant contigs in the HAP signature 6 days before the HAP onset.
#Input table : ALR/Presence absence counts

#In shell, run :
#format_input.py Contig_ALR_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("rpkm.res")

HAP_signature_contigs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5) +
  coord_flip() +
  theme_bw()  +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey")) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "LEfSe of contigs in HAP signature", x = "Differential abundant contigs", y = "LDA score (Log10)")+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text = element_text(size=30),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        plot.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = "bottom")


pdf("HAP_signature_contigs_LDA_PREVH.pdf",width=16,height=11);
HAP_signature_contigs
dev.off()

# FIGURE 5D - Phage lifestyle Barplot

#### Vibrant results

#VIBRANT to predict lifestyles (lysogenic/lytic) for contigs with minimum sequence length of 1000bp and containing at least 4 ORFs (open readings frames) 

singularity shell vibrant.sif
VIBRANT_run.py -i viral_contigs.fasta -t 114 -folder VIBRANT_results -virome

library(ggplot2)

# Define the data
data <- data.frame(
  Sample = c("IBIS", "PREVHAP", "IBIS", "PREVHAP"),
  Virus = c("Lytic", "Lytic", "Lysogenic", "Lysogenic"),
  Value = c(2582, 1385, 93, 67),
  Percent= c(96.3, 95.4, 3.7, 4.6)
)

# Define custom virus colors
style_colors <- c("Lytic" = "wheat1", "Lysogenic" = "orange3")

# Create the plot
plot <- ggplot(data, aes(fill = Virus, y = Percent, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = style_colors) +
  labs(x = "Cohort",
       y = "Percentage (%)") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  geom_text(aes(label = Value), position = position_stack(vjust = 0.5), size = 20, color = "black") +
  theme(legend.text = element_text(size = 40),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(face = "bold", size = 40)) +
  theme(axis.text.x = element_text(face = "bold", size = 40, colour = "black"))

# Print the plot
print(plot)


pdf("Vibrant.pdf",width=15,height=20);
plot
dev.off()


# FIGURE 5E - Caudoviricetes relative abundance in the sliding window 5-3 days before HAP onset

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("CAUDO_53_tracheo.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Caudoviricetes, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(Caudoviricetes ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Caudoviricetes ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

Caudoviricetes_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), Caudoviricetes)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) +  # Boxplot inside violin, adjust width here
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") +  # Exclude black from legend
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 100, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Caudoviricetes at 5-3 period
  before HAP onset (Relative abundance)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none")


pdf("Caudoviricetes_boxplot_53_tracheo.pdf",width=12,height=9);
Caudoviricetes_boxplot
dev.off()


# FIGURE 5F - Discriminant Caudoviricetes relative abundance in the sliding window 5-3 days before HAP onset


library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Caudo_53_lda_fish_tracheo.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Caudoviricetes, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(Caudoviricetes ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Caudoviricetes ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

Caudoviricetes_discriminant_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), Caudoviricetes)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) +  # Boxplot inside violin, adjust width here
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") +  # Exclude black from legend
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 25, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Discriminant Caudoviricetes at 5-3 period
  before HAP onset (Relative abundance)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none")


pdf("Caudoviricetes_discriminant_boxplot_53_tracheo.pdf",width=12,height=9);
Caudoviricetes_discriminant_boxplot
dev.off()

# FIGURE 5G Chordiagram HAP

#Run MAASLIN2 on RPKM counts of HAP-associated contigs with relative abundance of the core respiratory bacteriome

library(Maaslin2)

df_input_data = read.table(file = "RPKMcounts_HAP_associated.txt",
                          header = TRUE,
                          sep              = "\t", 
                          row.names        = 1,
                          stringsAsFactors = FALSE)

df_input_metadata = read.table(file = "16S_Core_respiratory_bacteriome.txt", 
                              header           = TRUE, 
                              sep              = "\t", 
                              row.names        = 1,
                              stringsAsFactors = FALSE)


fitData.ctrls.lefse= Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata,
  output = "Correlations_between_viral_contigs_and_core_bacteriome",
  normalization = "none",
  transform = "none",
  random_effects = c("Samples"))


#Use Maaslin2 output in Chordiagram

data <- read.delim("Correlations_viral_Core_Bacteriome_HAP.txt", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(circlize)
library(scales)
library(reshape)

grid.col = c(Fusobacterium = "red", Haemophilus = "orange", Prevotella = "green", Streptococcus = "blue", Veillonella = "purple")
col_fun = colorRamp2(c(-10,-5,0, 5, 10), c("red", "pink", "white", "lightblue","blue"))

pdf("chordDiagram_HAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()

# FIGURE 5H Chordiagram NO HAP

data <- read.delim("Correlations_viral_Core_Bacteriome_NOHAP.txt", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(circlize)
library(scales)
library(reshape)

grid.col = c(Fusobacterium = "red", Haemophilus = "orange", Prevotella = "green", Streptococcus = "blue", Veillonella = "purple")
col_fun = colorRamp2(c(0,0.2,0.4), c("white", "lightblue","blue"))

pdf("chordDiagram_noHAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()














##############################################################