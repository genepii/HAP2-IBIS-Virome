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

# FIGURE 5B - Betadiversity WBC/Hellinger/Sorensen in sliding window 5-3 days before HAP onset

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

##############################################################