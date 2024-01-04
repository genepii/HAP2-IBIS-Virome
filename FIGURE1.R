############################################################################################################################
#### CODE USED TO GENERATE FIGURE 1
############################################################################################################################

# ALR heatmap FIGURE 1A

library(pheatmap)
matrix <- read.delim("RPKMcounts.txt", header = TRUE, sep = "\t", row.names = 1)
matrix <- na.omit(matrix)
transformed_matrix <- log10(matrix+1)
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


# RELAB heatmap FIGURE 1B

library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)

matrix <- read.delim("Relative_abundance.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Family","Sample", "Value")

# Create a stacked bar plot with manual colors

virus_colors <- c("#450659FF","#460B5EFF","#471063FF","#481668FF","#481A6CFF","#481E70FF","#482374FF","#482778FF",
                  "#472C7AFF","#472F7EFF","#463480FF","#453882FF","#433D84FF","#424086FF","#0000FF","#3F4889FF","#3D4D8AFF",
                  "#3C508BFF","#3A538BFF","#38588CFF","#375B8DFF","#355F8DFF","#33628DFF","#31668EFF","#30698EFF","#2E6D8EFF",
                  "#2D708EFF","#2C738EFF","#2A768EFF","#29798EFF","#287D8EFF","#27808EFF","#25838EFF","#24868EFF","#238A8DFF",
                  "#228D8DFF","#21908CFF","#20938CFF","#1F968BFF","#1F9A8AFF","#1E9D89FF","#1FA188FF","#20A386FF","#22A785FF",
                  "#24AA83FF","#27AD81FF","#2AB07FFF","#2EB37CFF","#34B679FF","#39B977FF","#3FBC73FF","#44BF70FF","#4BC26CFF",
                  "#52C569FF","#59C864FF","#60CA60FF","#68CD5BFF","#70CF57FF","#78D152FF","#80D34DFF","#89D548FF","#92D742FF",
                  "#9BD93CFF","#A3DB36FF","#ADDC30FF","#B6DE2AFF","#C0DF25FF","#C9E020FF","#D2E21BFF","#DBE319FF", "#999999","#E4E419FF",
                  "#ECE51BFF","#F5E61FFF","#FDE725FF")


RELAB_heatmap <- ggplot(melted_data, aes(fill = Family, y = Value, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = virus_colors) +
  labs(title = "Relative abundance of viral families",
       x = "Samples",
       y = "Relative abundance") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_blank())

pdf("RELAB_heatmap.pdf",width=20,height=15);
RELAB_heatmap
dev.off()

# BETADIV 

# BRAY CURTIS

library(vegan)
data1<-read.delim("RPKMcounts.txt", row.names = 1)
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


***(the input file of betadiversity analysis contain the metric medians calculated between patients or within patients)***
	- Between analysis : we calculated the median of distance between each patient with others in HAP or no HAP groups
	- Within analysis : we calculated the median of distance between samples collected at different times and belonging to the same patient in HAP or no HAP groups 

"Tables containing median values from WBC, Hellinger or Sorensen served as input for statistical and dynamics analysis"
		
# FIGURE 1C - WBC (Weighted Bray-Curtis) HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_between_within_HAP.txt", stringsAsFactors = FALSE)
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

WBC_Between_within_violin_HAP <- ggplot(data, aes(GROUP, WBC)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) + 
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Within-HAP patients" = "orange"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Bray-Curtis dissimilarity"
  )

pdf("WBC_Between_within_violin_HAP.pdf",width=12,height=9);
WBC_Between_within_violin_HAP
dev.off()

# BETADIV FIGURE 1C - WBC no HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_between_within_NOHAP.txt", stringsAsFactors = FALSE)
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

WBC_Between_within_violin_noHAP <- ggplot(data, aes(GROUP, WBC)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) + 
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("Between-NO HAP patients" = "blue", "Within-NO HAP patients" = "lightblue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Bray-Curtis dissimilarity"
  )

pdf("WBC_Between_within_violin_noHAP.pdf",width=12,height=9);
WBC_Between_within_violin_noHAP
dev.off()

# BETADIV FIGURE 1C - Hellinger HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_between_within_HAP.txt", stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Within-HAP patients" = "orange"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.5, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Hellinger dissimilarity"
  )

pdf("Hellinger_Between_within_violin_HAP.pdf",width=12,height=9);
Hellinger_Between_within_violin_HAP
dev.off()

# BETADIV FIGURE 1C - Hellinger no HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_between_within_NOHAP.txt", stringsAsFactors = FALSE)
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

Hellinger_Between_within_violin_noHAP <- ggplot(data, aes(GROUP, Hellinger)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("Between-NO HAP patients" = "blue", "Within-NO HAP patients" = "lightblue"), name = "GROUP") +
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.5, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Hellinger dissimilarity"
  )

pdf("Hellinger_Between_within_violin_noHAP.pdf",width=12,height=9);
Hellinger_Between_within_violin_noHAP
dev.off()

# BETADIV FIGURE 1C - Sorensen HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_between_within_HAP.txt", stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Within-HAP patients" = "orange"), name = "GROUP") +
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Sorensen dissimilarity"
  )

pdf("Sorensen_Between_within_violin_HAP.pdf",width=12,height=9);
Sorensen_Between_within_violin_HAP
dev.off()

# BETADIV FIGURE 1C - Sorensen no HAP

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_between_within_NOHAP.txt", stringsAsFactors = FALSE)
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

Sorensen_Between_within_violin_noHAP <- ggplot(data, aes(GROUP, Sorensen)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) + 
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("Between-NO HAP patients" = "blue", "Within-NO HAP patients" = "lightblue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Sorensen dissimilarity"
  )

pdf("Sorensen_Between_within_violin_noHAP.pdf",width=12,height=9);
Sorensen_Between_within_violin_noHAP
dev.off()