############################################################################################################################
#### CODE USED TO GENERATE FIGURE 1
############################################################################################################################

# ALR heatmap FIGURE 1A

library(pheatmap)
matrix <- read.delim("counts.txt", header = TRUE, sep = "\t", row.names = 1)
matrix <- na.omit(matrix)
transformed_matrix <- log(matrix+0.01)
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
pdf("ALR_heatmap.pdf",width=20,height=14);
ALR_heatmap
dev.off()

# FIGURE 1B - Phages lifestyles

#VIBRANT to predict lifestyles (lysogenic/lytic) for contigs with minimum sequence length of 1000bp and containing at least 4 ORFs (open readings frames) 

singularity shell vibrant.sif
VIBRANT_run.py -i viral_contigs.fasta -t 114 -folder VIBRANT_results -virome

library(ggplot2)

# Define the data
data <- data.frame(
  Sample = c("HAP", "NO HAP", "HAP", "NO HAP"),
  Virus = c("Lytic", "Lytic", "Lysogenic", "Lysogenic"),
  Value = c(1734, 1648, 74, 72),
  Percent= c(95.9, 95.8, 4.1, 4.2)
)

# Define custom virus colors
style_colors <- c("Lytic" = "wheat1", "Lysogenic" = "orange3")

# Create the plot
plot <- ggplot(data, aes(fill = Virus, y = Percent, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = style_colors) +
  labs(x = "",
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

pdf("Vibrant_IBIS.pdf",width=15,height=20);
plot
dev.off()

# FIGURE 1C - BETADIVERSITY

# BETADIV 

# BRAY CURTIS

library(vegan)
data1<-read.delim("counts.txt", row.names = 1)
dataTransposed1<-t(data1)
dis <- vegdist(dataTransposed1, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"WeightedBrayCurtis.txt", sep = '\t')

# Hellinger

library(adespatial)
data2 <-read.delim("ALRcounts.txt", row.names=1)
dataTransposed2 <-t(data2)
dist.hel <-dist.ldc(dataTransposed2, method = "hellinger")
dist2 <-as.matrix(dist.hel)
write.table(dist2, "Hellinger.txt", sep = '\t')

# Sorensen 

library(vegan)
data3<-read.delim("Presence_absence_counts.txt", row.names = 1)
dataTransposed3<-t(data3)
dis <- vegdist(dataTransposed3, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"Sorensen.txt", sep = '\t')

# Reformat Betadiv superimposed Output table (BRAY CURTIS, HELLINGER and Sorensen)

library(reshape2)

data <- read.delim("Betadiv_superimposed_matrix.txt", stringsAsFactors = FALSE, row.names = 1, header = TRUE, check.names = FALSE)

superimposed_matrix <- as.matrix(data)
upper_logical <- upper.tri(superimposed_matrix)
upper_triangle_subset <- matrix(NA, nrow = nrow(superimposed_matrix), ncol = ncol(superimposed_matrix))
upper_triangle_subset[upper_logical] <- superimposed_matrix[upper_logical]
row_names <- row.names(data)
col_names <- colnames(data)
rownames(upper_triangle_subset) <- row_names
colnames(upper_triangle_subset) <- col_names
print(upper_triangle_subset)

# Reshape the data into a vertical table format
melted_data <- melt(upper_triangle_subset, id.vars = "sample1", variable.name = "sample2", value.name = "value")
colnames(melted_data) <- c("sample1", "sample2", "value")
write.table(melted_data, file="Formated_betadiv_table.txt", sep = "\t")


***(the input file of betadiversity analysis contain the metric medians calculated between patients or within patients)***
	- Inter/Between analysis : we calculated the median of distance between each patient with others
	- Intra/Within analysis : we calculated the median of distance between samples collected at different times and belonging to the same patient 

"Tables containing median values from WBC, Hellinger or Sorensen served as input for statistical and dynamics analysis"

# FIGURE 1C - WBC

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_between_within.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(WBC, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(WBC ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(WBC ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

WBC_Between_within_violin <- ggplot(data, aes(GROUP, WBC)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-patients" = "magenta", "Within-patients" = "green3"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Bray-Curtis dissimilarity"
  )

pdf("WBC_Between_within_violin.pdf",width=12,height=9);
WBC_Between_within_violin
dev.off()


# BETADIV FIGURE 1C  - Hellinger

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_between_within.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Hellinger, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(Hellinger ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Hellinger ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

Hellinger_Between_within_violin_HAP <- ggplot(data, aes(GROUP, Hellinger)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-patients" = "magenta", "Within-patients" = "green3"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.5, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Hellinger dissimilarity"
  )

pdf("Hellinger_Between_within_violin.pdf",width=12,height=9);
Hellinger_Between_within_violin_HAP
dev.off()



# BETADIV FIGURE 1C  - Sorensen

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_between_within.txt", stringsAsFactors = FALSE)
head(data)

data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Sorensen, type = "median_iqr")

stat.test <- data %>% 
  wilcox_test(Sorensen ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Sorensen ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

Sorensen_Between_within_violin_HAP <- ggplot(data, aes(GROUP, Sorensen)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-patients" = "magenta", "Within-patients" = "green3"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Sorensen dissimilarity"
  )

pdf("Sorensen_Between_within_violin.pdf",width=12,height=9);
Sorensen_Between_within_violin_HAP
dev.off()
	
