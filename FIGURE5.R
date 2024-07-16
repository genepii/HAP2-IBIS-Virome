############################################################################################################################
#### CODE USED TO GENERATE FIGURE 5
############################################################################################################################

# FIGURE 5A

library(gplots)
library(ComplexHeatmap)

# Read the matrix and keep NA values
matrix <- read.delim("RPKMcounts.txt", header = TRUE, sep = "\t", row.names = 1)

# Perform log transformation, handling zero values to avoid -Inf
transformed_matrix <- log10(matrix)
transformed_matrix[transformed_matrix=="-Inf"] <- NA

# Read the sample annotations
my_sample_col <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)

# Select top 10 families based on row sums
family_sums <- rowSums(transformed_matrix, na.rm = TRUE)
top10_families <- names(sort(family_sums, decreasing = TRUE))[1:10]
matrix_top10 <- transformed_matrix[top10_families, ]


col = list(HAP.condition = c("HAP" = "red", "no HAP" = "blue"),
           ARDS.condition = c("ARDS" = "green", "no ARDS" = "purple"),
           Death = c("Deceased" = "yellow", "Alive" = "orange"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  HAP.condition = my_sample_col$HAP.condition, 
  ARDS.condition = my_sample_col$ARDS.condition, 
  Death = my_sample_col$Death, 
  HAP.onset = my_sample_col$HAP.onset,
  col = col,
  annotation_legend_param = list(
    HAP.condition = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 30)),
    ARDS.condition = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 30)),
    Death = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 30))
  ),
  annotation_name_gp = gpar(fontsize = 30)
)


pdf("heatmap.pdf", width = 20, height = 14)

# Combine the heatmap and the annotation
Heatmap(matrix_top10, name = "log10RPKM",
        top_annotation = ha,
        row_names_gp = gpar(fontsize = 30),
        column_names_gp = gpar(fontsize = 0))

dev.off()

# FIGURE 5B - Phage lifestyle Barplot

#### Vibrant results

#VIBRANT to predict lifestyles (lysogenic/lytic) for vOTUs with minimum sequence length of 1000bp and containing at least 4 ORFs (open readings frames) 

singularity shell vibrant.sif
VIBRANT_run.py -i viral_vOTUs.fasta -t 114 -folder VIBRANT_results -virome

library(ggplot2)

# Define the data (Vibrant output)

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

# FIGURE 5C - Betadiversity WBC/Hellinger/Sorensen in sliding window 5-3 days before HAP onset

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Betadiv_metric.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Betadiv, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(Betadiv ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Betadiv ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

Betadiv_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), Betadiv)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 100, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Betadiv") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none")


pdf("Betadiv_boxplot_53.pdf",width=12,height=9);
Betadiv_boxplot
dev.off()

# FIGURE 5D

#Use Fisher test to identify discriminant vOTUs in HAP and no HAP signature 5-3 days before the HAP onset.
#Input table : Presence absence counts

otu_data <- read.delim("Presence absence.txt", header = TRUE)

for (n in 1:nrow(otu_data)){
  otu_data_2<-otu_data[n,-1]
  otu_data_3<-rbind(otu_data_2,apply(otu_data[-n,-1],2,FUN=sum))
  otu_data$p.value[n]<-fisher.test(otu_data_3)$p.value
}

library(ggplot2)

data <- read.delim("fisher_output.txt")

HAP_noHAP_signature_vOTUs <- ggplot(data, aes(reorder(OTU, PVAL), PVAL, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.8, size = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey28", "Eukaryotic viruses" = "grey")) +
  theme(strip.placement = "outside",
                     strip.text.y = element_text(angle = 0)) +
  labs(title = "", x = "Significant vOTUs", y = "Fisher test -log10(P-Value)")+
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


pdf("signature_Fisher_PREVHAP.pdf",width=16,height=11);
HAP_noHAP_signature_vOTUs
dev.off()


# FIGURE 5E - Lefse HAP vOTUs

#Use LEfSe to identify discriminant vOTUs in the HAP signature 5-3 days before the HAP onset.
#Input table : log10RPKM counts

#In shell, run :
#format_input.py vOTUs_log10RPKM_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("LDA.res")

HAP_signature_vOTUs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5) +
  coord_flip() +
  theme_bw()  +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey28", "Eukaryotic viruses" = "grey")) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "LEfSe of vOTUs in HAP signature", x = "Differential abundant vOTUs", y = "LDA score (Log10)")+
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


pdf("HAP_signature_vOTUs_LDA_PREVH.pdf",width=16,height=11);
HAP_signature_vOTUs
dev.off()




# FIGURE 5F-G 

#Run correlation tests on  log10(relative_abundace) of HAP-associated vOTUs and the core respiratory bacteriome

# Load necessary libraries
library(corrplot)
library(dplyr)
library(Hmisc)
library(tibble)

# Read the data from text files
table1 <- read.table("log10(relab)_16S.txt", header = TRUE)
table2 <- read.table("log10(relab)_vOTUs.txt", header = TRUE)

# Transpose table2 and convert row names to a column
table2 <- as.data.frame(t(table2))
table2 <- rownames_to_column(table2, var = "Samples")

# Rename the first column to match the merging key in table1
colnames(table2)[1] <- "Samples"

# Merge the data frames by Sample ID
merged_data <- inner_join(table1, table2, by = "Samples")

# Remove the SampleID column for correlation calculation
merged_data <- merged_data %>% select(-Samples)

# Compute the Spearman correlation matrix and p-values
cor_res <- rcorr(as.matrix(merged_data), type = "spearman")
cor_matrix <- cor_res$r
p_matrix <- cor_res$P

# Export the correlation matrix and p-value matrix to text files
write.table(cor_matrix, file = "correlation_matrix.txt", sep = "\t", quote = FALSE)
write.table(p_matrix, file = "pvalue_matrix.txt", sep = "\t", quote = FALSE)


#Use correlation output in Chordiagram

data <- read.delim("correlations.txt", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(circlize)
library(scales)
library(reshape)

grid.col = c(Fusobacterium = "red", Haemophilus = "orange", Prevotella = "green", Streptococcus = "blue", Veillonella = "purple")
col_fun = colorRamp2(c(-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7,0.8,0.9,1), c("red4","darkred","coral4","red3","red2","red1","red","slateblue", "slateblue1", "slateblue2", "slateblue3", "slateblue4","darkblue","navy"))

pdf("chordDiagram_noHAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid"), col = col_fun, grid.col = grid.col, scale = F)
lgd <- Legend(labels = c(-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,0.4,0.5,0.6,0.7,0.8,0.9,1), 
              legend_gp = gpar(fill = c("red4","darkred","coral4","red3","red2","red1","red","slateblue", "slateblue1", "slateblue2", "slateblue3", "slateblue4","darkblue","navy")), 
              title = "Correlation",
              title_position = "leftcenter-rot",
              labels_rot = TRUE,
              grid_height = unit(1, "cm"),  
              grid_width = unit(1, "cm"),   
              title_gp = gpar(fontsize = 50),  
              labels_gp = gpar(fontsize = 35))
draw(lgd, x = unit(1, "npc") - unit(1, "cm"), y = unit(1, "npc") - unit(1, "cm"), just = c("right", "top"))
dev.off()

##############################################################