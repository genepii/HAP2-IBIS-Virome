############################################################################################################################
#### CODE USED TO GENERATE FIGURE 4
############################################################################################################################

# Statistical test for dynamics analysis: Mann-Whitney U Test and bootstrapping

#betadiv_metric can be a Bray-Curtis, Hellinger or Sorensen median values

library(boot)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

dat <- read.delim("betadiv_medians.txt", stringsAsFactors = FALSE)
result<-wilcox.test(betadiv_metric ~ GROUP, data=dat, alternative="two.sided", correct=TRUE,conf.int=TRUE, conf.level=0.95)
print(result)

Mboot = boot(dat$betadiv_metric,
             function(x,i) median(x[i]),
             R=10000)
boot.ci(Mboot,
        conf = 0.95,
        type = "norm")

# FIGURE 4A - Dynamics WBC

data <- read.delim("WBC_sliding_window.txt", stringsAsFactors = FALSE)

data$GROUP <- factor(data$GROUP, levels = c("upcoming HAP", "no HAP"))

custom_colors <- c("upcoming HAP" = "pink", "no HAP" = "blue")

WBC_dynamics <- ggplot(data, aes(x = DAY, y = WBC, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = ifelse(data$GROUP == "no HAP", 0.3, 0.6)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size=25),
    legend.title = element_blank(),
    legend.position = "none"
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Bray Curtis Dissimilarity")


WBC_dynamics_stat <- WBC_dynamics + scale_y_continuous(breaks = c(0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("5-4","4-3", "3-2", "2-1", "1-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -16.8)

pdf("WBC_dynamics.pdf",width=15,height=10);
WBC_dynamics_stat
dev.off()




# FIGURE 4A - Dynamics Hellinger

data <- read.delim("Hellinger_sliding_window.txt", stringsAsFactors = FALSE)

data$GROUP <- factor(data$GROUP, levels = c("upcoming HAP", "no HAP"))

custom_colors <- c("upcoming HAP" = "pink", "no HAP" = "blue")

Hellinger_dynamics <- ggplot(data, aes(x = DAY, y = Hellinger, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = ifelse(data$GROUP == "no HAP", 0.3, 0.6)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) + 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size=25),
    legend.title = element_blank(),
    legend.position = "none"
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Hellinger Dissimilarity")


Hellinger_dynamics_stat <- Hellinger_dynamics + scale_y_continuous(breaks = c(0.8,0.9,1,1.1,1.2,1.3,1.4)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("5-4","4-3", "3-2", "2-1", "1-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -16.5)

pdf("Hellinger_dynamics.pdf",width=15,height=10);
Hellinger_dynamics_stat
dev.off()


# FIGURE 4A - Dynamics Sorensen

data <- read.delim("Sorensen_sliding_window.txt", stringsAsFactors = FALSE)

data$GROUP <- factor(data$GROUP, levels = c("upcoming HAP", "no HAP"))

custom_colors <- c("upcoming HAP" = "pink", "no HAP" = "blue")

Sorensen_dynamics <- ggplot(data, aes(x = DAY, y = Sorensen, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = ifelse(data$GROUP == "no HAP", 0.3, 0.6)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) +  
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size=25),
    legend.title = element_blank(),
    legend.position = "none"
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Sorensen Dissimilarity")


Sorensen_dynamics_stat <- Sorensen_dynamics + scale_y_continuous(breaks = c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("5-4","4-3", "3-2", "2-1", "1-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -6.5)

pdf("Sorensen_dynamics.pdf",width=15,height=10);
Sorensen_dynamics_stat
dev.off()

# FIGURE 4B - Fisher HAP/no HAP vOTUs

#Use Fisher test to identify discriminant vOTUs signature 5-4 days before the HAP onset.
#Input table : Presence absence counts

otu_data <- read.delim("presence_absence.txt", header = TRUE)

for (n in 1:nrow(otu_data)){
  otu_data_2<-otu_data[n,-1]
  otu_data_3<-rbind(otu_data_2,apply(otu_data[-n,-1],2,FUN=sum))
  otu_data$p.value[n]<-fisher.test(otu_data_3)$p.value
}

# FIGURE 4C - Lefse

#Use LEfSe to identify discriminant vOTUs in the signature 5-4 days before the HAP onset.
#Input table : log10RPKM counts

#In shell, run :
#format_input.py vOTUs_log10RPKM_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res

#Merge Lefse and fisher outputs in one table and plot heatmap

library(pheatmap)
matrix <- read.delim("data.txt", header = TRUE, sep = "\t", row.names = 1)
my_sample_col <- read.table("metacol.txt", header = TRUE, sep = "\t",row.names = 1)
my_sample_row <- read.table("metarow.txt", header = TRUE, sep = "\t",row.names = 1)
ann_colors = list(
  Group = c("upcoming HAP" = "pink", "no HAP" = "blue")
)

heatmap <- pheatmap(
  matrix,
  scale = "row",
  color = colorRampPalette(c('darkblue', 'white', 'red'))(100),
  annotation_row = my_sample_row,
  annotation_col = my_sample_COL,
  cluster_cols = F,
  cellheight = 10,
  cellwidth = 5,
  cluster_rows =T,
  annotation_colors = ann_colors,
  fontsize = 10, 
  fontsize_row = 10,
  fontsize_col = 5, 
  legend = TRUE,
  border_color = "white",
  na_col = "grey",
  annotation_legend = TRUE, show_colnames = F, show_rownames = F, angle_col = 90
)
pdf("heatmap.pdf",width=20,height=14)
heatmap
dev.off()

# FIGURE 4D

#Run correlation tests on  log10(relative_abundace) of viral signature vOTUs and the core respiratory bacteriome

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


#Use correlation output for bubble plot

library(ggtext)
library(dplyr)
library(ggplot2)

results_sig = read.delim(file="Correlations.txt",header=T,check.names=F)

results_sig <- results_sig %>%
  arrange(desc(Occurrence_Percentage))  # Sort the data

results_sig$OTU <- factor(results_sig$OTU, levels = unique(results_sig$OTU))

correlations <- ggplot(results_sig, aes(x = row, y = OTU)) +
  geom_point(aes(size = abs(Occurrence_Percentage), color = cor, alpha = 1, shape = Lifestyle)) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue") +
  scale_shape_manual(values = c("lytic" = 16, "lysogenic" = 17, "Unknown" = 18)) +
  coord_flip() +
  labs(
    x = "", 
    y = "Feature", 
    size = "Occurrence Percentage"
  ) +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 1, angle = 90, hjust = 1),
    axis.text.y = element_markdown(size = 12)
  ) + facet_grid(Type~group, scales = "free")

correlations

pdf("correlations.pdf",width=30,height=10);
correlations
dev.off()
