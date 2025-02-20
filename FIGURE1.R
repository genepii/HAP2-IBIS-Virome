############################################################################################################################
#### CODE USED TO GENERATE FIGURE 1
############################################################################################################################
# Effect size 1A-B

library(vegan)
library(ggplot2)

# Read data
otu = read.delim(file="dRPKM.txt", header=T, row.names=1, check.names=F)
env = read.delim(file="metadata.txt", header=T, row.names = 1, check.names=F) 

# Hellinger transformation of the OTU data
otu_hel <- decostand(otu, method = 'hellinger')
otu_hel = t(otu_hel)

# Perform PCA
otupca <- rda(otu_hel, env, scale = FALSE)
summary(otupca, scaling = 1)
plot(otupca, choices = c(1, 2), scaling = 1,  display = c('wa'))

# Perform envfit to assess the relationship between environmental variables and the PCA axes
otupca_ef <- envfit(otupca, env, perm = 999, choices = c(1, 2), display = 'sites')
otupca_ef

# p-value adjustment with FDR method
otupca_ef_adj <- otupca_ef
otupca_ef_adj$factors$pvals <- p.adjust(otupca_ef_adj$factors$pvals, method = 'fdr')
otupca_ef_adj$vectors$pvals <- p.adjust(otupca_ef_adj$vectors$pvals, method = 'fdr')
otupca_ef_adj

env_effect_size <- otupca_ef$vectors$r
env_variables <- colnames(env)
env_p_values <- otupca_ef$vectors$pvals 

effect_size_df <- data.frame(
  Variable = env_variables,
  EffectSize = env_effect_size,
  PValue = env_p_values
)

effect_size_df$PercentageVariation <- (effect_size_df$EffectSize / sum(effect_size_df$EffectSize)) * 100

effect_size_df$AdjPValue <- p.adjust(effect_size_df$PValue, method = "fdr")

significance_threshold <- 0.05

effect_size_df$Significant <- ifelse(effect_size_df$PValue < significance_threshold, 
                                     "*", "")

effect_size_df$Significant <- ifelse(effect_size_df$PValue < 0.001, "***", effect_size_df$Significant)
effect_size_df$Significant <- ifelse(effect_size_df$PValue < 0.01 & effect_size_df$PValue >= 0.001, "**", effect_size_df$Significant)
effect_size_df$Significant <- ifelse(effect_size_df$PValue < 0.05 & effect_size_df$PValue >= 0.01, "*", effect_size_df$Significant)

# Plot the effect sizes with significance stars using ggplot2
ES<-ggplot(effect_size_df, aes(x = reorder(Variable, EffectSize), y = EffectSize, fill = Variable)) +
  geom_bar(stat = 'identity') +  # Plot the effect sizes
  geom_text(aes(label = Significant), vjust = -0.5, hjust = 0, size = 10) +
  coord_flip() +  
  xlab('') +
  ylab('Effect Size (R²)') +
  theme_classic() +
  scale_fill_brewer(palette = "Set3")+
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        legend.position = "none")

ES
pdf("effect_size.pdf", height = 7, width = 12)
ES
dev.off()

# heatmap 1C

library(gplots)
library(ComplexHeatmap)

# Read the matrix and keep NA values
matrix <- read.delim("dRPKM.txt", header = TRUE, sep = "\t", row.names = 1)

# Read the sample annotations
my_sample_col <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)

matrix<-matrix[apply(matrix,1,function(x)length(x[x>0]))>=0.5*ncol(matrix),]

matrix <- matrix[, colSums(matrix) != 0]

# Keep only the samples present in the matrix
my_sample_col <- my_sample_col[colnames(matrix), , drop = FALSE]

# Select top 30 genera based on row sums
family_sums <- rowSums(matrix, na.rm = TRUE)
top10_families <- names(sort(family_sums, decreasing = TRUE))[1:30]
matrix_top10 <- matrix[top10_families, ]

# Perform log transformation, handling zero values to avoid -Inf
transformed_matrix <- log10(matrix)
transformed_matrix[transformed_matrix=="-Inf"] <- NA

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
  annotation_name_gp = gpar(fontsize = 10)
)


pdf("heatmap.pdf", width = 20, height = 14)
 
Heatmap(
  transformed_matrix, 
  name = "log10RPKM",
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha,  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 0),
  heatmap_legend_param = list(
    legend_direction = "horizontal", 
    legend_position = "top"          
  )
)

dev.off()


##### Otu prevalence_abundance 1D

library(ggplot2)
library(reshape2)
library(readr)
library(RColorBrewer)

data  = read.delim(file="data.txt",header=T,check.names=F)

data$Occurrence_Percentage <- as.numeric(data$Occurrence_Percentage)

sig <- ggplot(data, aes(x = MRA, y = Occurrence_Percentage, color = Type)) +
  geom_point(alpha = 0.4, size = 5) +
  theme_classic() +
  labs(title = "", x = "Log10MRA", y = "Occurrence Percentage") +
  theme(
    legend.position = "right"
    ,
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
  )

sig

pdf("MRA_occurence.pdf",width=5,height=5);
sig
dev.off()

#1E - Phages lifestyles
in bash run phabox2
#conda activate phabox2
#phabox2 --task end_to_end --dbdir /phabox_db_v2 --outpth  output_folder --contigs vOTUs.fasta --threads 70 --len 300

library(ggplot2)

library(dplyr)

# Define the raw dataset
dara_raw=phabox_predictions
# Convert to long format for ggplot
data_long <- tidyr::pivot_longer(data_raw, cols = c("Lytic", "Lysogenic"), names_to = "Virus", values_to = "Value")

# Calculate percentages
data_long <- data_long %>%
  group_by(Sample) %>%
  mutate(Percent = (Value / sum(Value)) * 100)

# Define custom virus colors
style_colors <- c("Lytic" = "gold2", "Lysogenic" = "seagreen")

# Create the plot
plot <- ggplot(data_long, aes(fill = Virus, y = Percent, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = style_colors) +
  labs(x = "",
       y = "Percentage (%)") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  geom_text(aes(label = Value), position = position_stack(vjust = 0.5), size = 6, color = "black") +
  theme(legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14, colour = "black"))

# Print the plot
print(plot)

pdf("phabox.pdf",width=7.5,height=10)
plot
dev.off()

# FIGURE 1F - BETADIVERSITY

# BETADIV 

# BRAY CURTIS

library(vegan)
data1<-read.delim("RPKM_bacteriome_deconta_rarefied.txt", row.names = 1)
dataTransposed1<-t(data1)
dis <- vegdist(dataTransposed1, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"WeightedBrayCurtis.txt", sep = '\t')

# Hellinger

library(adespatial)
data2 <-read.delim("RPKM_bacteriome_deconta_rarefied.txt", row.names=1)
data2 <-log10(data2+1)
dataTransposed2 <-t(data2)
dist.hel <-dist.ldc(dataTransposed2, method = "hellinger")
dist2 <-as.matrix(dist.hel)
write.table(dist2, "Hellinger.txt", sep = '\t')

# Sorensen 

library(vegan)
data3<-read.delim("RPKM_bacteriome_deconta_rarefied.txt", row.names = 1)
data3[data3 > 0] <- 1
dataTransposed3<-t(data3)
dis <- vegdist(dataTransposed3, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"Sorensen.txt", sep = '\t')

# Reformat Betadiv superimposed Output table (BRAY CURTIS, HELLINGER and Sorensen)

library(reshape2)

data <- read.delim("Hellinger.txt", stringsAsFactors = FALSE, row.names = 1, header = TRUE, check.names = FALSE)

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
write.table(melted_data, file="hel_format.txt", sep = "\t")


***(the input file of betadiversity analysis contain the metric medians calculated between patients or within patients)***
	- Inter/Between analysis : we calculated the median of distance between each patient with others
	- Intra/Within analysis : we calculated the median of distance between samples collected at different times and belonging to the same patient 

"Tables containing median values from WBC, Hellinger or Sorensen served as input for statistical and dynamics analysis"

# FIGURE 1F - WBC

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

dat <- read.delim("WBC.txt", stringsAsFactors = FALSE)
head(dat)

dat %>% sample_n_by(GROUP, size = 2)
dat %>%
  group_by(GROUP) %>%
  get_summary_stats(VirF, type = "median_iqr")

stat.test <- dat %>% 
  wilcox_test(VirF ~ GROUP) %>%
  add_significance()
stat.test
dat %>% wilcox_effsize(VirF ~ GROUP)
stat.test <- stat.test %>% add_xy_position(x = "GROUP")

ggplot(dat, aes(GROUP, VirF)) +  
  # Boxplot layer: will be drawn first, in the background
  geom_boxplot(aes(fill = GROUP), width = 0.4, color = "black", outlier.shape = NA, alpha = 0.2) +  
  scale_fill_manual(values = c("C1" = "red", "C2" = "blue"), name = "GROUP") +  
  
  # Add p-value annotation with stat_pvalue_manual()
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 9, bracket.size = 2) +
  
  # Customize labels and theme
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40),
    axis.text.x = element_text(size = 40, colour = "black"),
    axis.text.y = element_text(size = 40),
    plot.subtitle = element_text(size = 30),
    legend.text = element_text(size = 40),
    legend.position = "none"
  ) +
  
  # Y-axis label
  labs(y = "Shannon_BAC") +
  
  # Jittered points layer: will be drawn last, in the foreground
  geom_jitter(aes(color = GROUP), width = 0.03, size = 10, alpha = 0.4) +
  scale_color_manual(values = c("C1" = "red", "C2" = "blue"))


#pdf("WBC_Between_within_violin.pdf",width=12,height=9);
#WBC_Between_within_violin
#dev.off()


# BETADIV FIGURE 1F  - Hellinger

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
  scale_fill_manual(values = c("Inter patients" = "magenta", "Intra patients" = "green3"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.5, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Hellinger dissimilarity"
  )

pdf("Hellinger_Between_within_violin.pdf",width=12,height=9);
Hellinger_Between_within_violin_HAP
dev.off()



# BETADIV FIGURE 1F  - Sorensen

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
  scale_fill_manual(values = c("Inter patients" = "magenta", "Intra patients" = "green3"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Sorensen dissimilarity"
  )

pdf("Sorensen_Between_within_violin.pdf",width=12,height=9);
Sorensen_Between_within_violin_HAP
dev.off()
	
