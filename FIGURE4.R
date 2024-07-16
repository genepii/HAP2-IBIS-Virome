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

# FIGURE 4B - Caudoviricetes dynamics


library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Caudoviricetes_sliding_windows.txt", stringsAsFactors = FALSE)

data$GROUP <- factor(data$GROUP, levels = c("upcoming HAP", "no HAP"))

custom_colors <- c("upcoming HAP" = "pink", "no HAP" = "blue")

Caudoviricetes_dynamics <- ggplot(data, aes(x = DAY, y = Caudoviricetes, group = GROUP, color = GROUP)) +
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
    axis.title.y = element_text(size = 35),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    legend.text = element_text(size=25),
    legend.title = element_blank(),
    legend.position = "none"
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Caudoviricetes relative abundance (%)")


Caudoviricetes_dynamics_stat <- Caudoviricetes_dynamics + scale_y_continuous(breaks = c(60,70,80,90,100)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("5-4","4-3", "3-2", "2-1", "1-0"))+
  geom_text(data = data, aes(x = DAY, y = 60, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -14.3)

pdf("Caudoviricetes_dynamics.pdf",width=15,height=10);
Caudoviricetes_dynamics_stat
dev.off()


# FIGURE 4C - Fisher HAP/no HAP vOTUs

#Use Fisher test to identify discriminant vOTUs in HAP and no HAP signature 5-4 days before the HAP onset.
#Input table : Presence absence counts

otu_data <- read.delim("presence_absence.txt", header = TRUE)

for (n in 1:nrow(otu_data)){
  otu_data_2<-otu_data[n,-1]
  otu_data_3<-rbind(otu_data_2,apply(otu_data[-n,-1],2,FUN=sum))
  otu_data$p.value[n]<-fisher.test(otu_data_3)$p.value
}


library(ggplot2)

data <- read.delim("fisher_output.txt")

HAP_noHAP_specific_vOTUs <- ggplot(data, aes(reorder(OTU, PVAL), PVAL, fill = CLASS)) +
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
        strip.text = element_text(size=40),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = "bottom") + facet_wrap(~GROUP)


pdf("signature_vOTUs_fisher.pdf",width=16,height=11);
HAP_noHAP_specific_vOTUs
dev.off()

# FIGURE 4D - Lefse HAP vOTUs

#Use LEfSe to identify discriminant vOTUs in the HAP signature 5-4 days before the HAP onset.
#Input table : log10RPKM counts

#In shell, run :
#format_input.py vOTUs_log10RPKM_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("LDA_hapnohap.res")

HAP_discriminant_vOTUs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Other bacteriophages" ="#ffcc00","Caudoviricetes" = "pink2")) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "LEfSe of vOTUs in HAP signature", x = "Differential abundant vOTUs", y = "LDA score (Log10)")+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text = element_text(size=40),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        plot.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = "bottom") + facet_wrap(~Group)


pdf("HAP_lefse_output.pdf",width=16,height=11);
HAP_discriminant_vOTUs
dev.off()

# FIGURE 4E-F Chordiagram 

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

