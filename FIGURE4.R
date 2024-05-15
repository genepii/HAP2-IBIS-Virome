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
wilcox.test(dat$betadiv_metric,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)
Mboot = boot(dat$betadiv_metric,
             function(x,i) median(x[i]),
             R=10000)
boot.ci(Mboot,
        conf = 0.95,
        type = "norm")

# FIGURE 4A - Dynamics WBC

data <- read.delim("WBC_dynamics.txt", stringsAsFactors = FALSE)

custom_colors <- c("Upcoming HAP" = "pink", "NO_HAP" = "blue")

WBC_dynamics <- ggplot(data, aes(x = DAY, y = WBC, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = 0.4) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size=25),
    legend.title = element_blank()
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Bray Curtis Dissimilarity")


WBC_dynamics_stat <- WBC_dynamics + scale_y_continuous(breaks = c(0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("6-4","5-3", "4-2", "3-1", "2-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -16)

pdf("WBC_dynamics.pdf",width=15,height=10);
WBC_dynamics_stat
dev.off()


# FIGURE 4A - Dynamics Hellinger

data <- read.delim("Hellinger_dynamics.txt", stringsAsFactors = FALSE)

custom_colors <- c("Upcoming HAP" = "pink", "NO_HAP" = "blue")

Hellinger_dynamics <- ggplot(data, aes(x = DAY, y = Hellinger, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = 0.4) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) + 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size=25),
    legend.title = element_blank()
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Hellinger Dissimilarity")


Hellinger_dynamics_stat <- Hellinger_dynamics + scale_y_continuous(breaks = c(0.8,0.9,1,1.1,1.2,1.3,1.4)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("6-4","5-3", "4-2", "3-1", "2-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -17.5)

pdf("Hellinger_dynamics.pdf",width=15,height=10);
Hellinger_dynamics_stat
dev.off()


# FIGURE 4A - Dynamics Sorensen

data <- read.delim("Sorensen_dynamics.txt", stringsAsFactors = FALSE)

custom_colors <- c("Upcoming HAP" = "pink", "NO_HAP" = "blue")

Sorensen_dynamics <- ggplot(data, aes(x = DAY, y = Sorensen, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = 0.4) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) +  
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size=25),
    legend.title = element_blank()
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Sorensen Dissimilarity")


Sorensen_dynamics_stat <- Sorensen_dynamics + scale_y_continuous(breaks = c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("6-4","5-3", "4-2", "3-1", "2-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.8, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -6.8)

pdf("Sorensen_dynamics.pdf",width=15,height=10);
Sorensen_dynamics_stat
dev.off()

# FIGURE 4B - Dynamics Caudoviricetes

data <- read.delim("Caudoviricetes_dynamics.txt", stringsAsFactors = FALSE)

custom_colors <- c("Upcoming HAP" = "pink", "NO_HAP" = "blue")

Caudoviricetes_dynamics <- ggplot(data, aes(x = DAY, y = Caudoviricetes, group = GROUP, color = GROUP)) +
  geom_line() +
  geom_point(shape=20, size=10) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = GROUP), alpha = 0.4) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(-4,3)) +  
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size=25),
    legend.title = element_blank()
  ) +  
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Caudoviricetes relative abundance")


Caudoviricetes_dynamics_stat <- Caudoviricetes_dynamics + scale_y_continuous(breaks = c(60,70,80,90,100)) +scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("6-4","5-3", "4-2", "3-1", "2-0"))+
  geom_text(data = data, aes(x = DAY, y = 60, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -12.5)

pdf("Caudoviricetes_dynamics.pdf",width=15,height=10);
Caudoviricetes_dynamics_stat
dev.off()

# FIGURE 4C - Fisher HAP/no HAP contigs

#Use Fisher test to identify discriminant contigs in HAP and no HAP signature 6 days before the HAP onset.
#Input table : Presence absence counts

otu_data <- read.delim("otu.txt", header = TRUE)

for (n in 1:nrow(otu_data)){
  otu_data_2<-otu_data[n,-1]
  otu_data_3<-rbind(otu_data_2,apply(otu_data[-n,-1],2,FUN=sum))
  otu_data$p.value[n]<-fisher.test(otu_data_3)$p.value
}


library(ggplot2)

data <- read.delim("HAP_noHAP_fisher_output.txt")

HAP_noHAP_specific_contigs <- ggplot(data, aes(reorder(OTU, PVAL), PVAL, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.8, size = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey", "Eukaryotic viruses" = "red")) +
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


pdf("HAP_noHAP_signature_contigs_fisher.pdf",width=16,height=11);
HAP_noHAP_specific_contigs
dev.off()

# FIGURE 4D - Lefse HAP contigs

#Use LEfSe to identify discriminant contigs in the HAP signature 6 days before the HAP onset.
#Input table : ALR/Presence absence counts

#In shell, run :
#format_input.py Contig_ALR_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("LDA_hapnohap.res")

HAP_discriminant_contigs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = CLASS)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Caudoviricetes" = "pink2","Other bacteriophages" ="#ffcc00", "Unclassified viruses" = "grey")) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "LEfSe of contigs in HAP signature", x = "Differentially abundant contigs", y = "LDA score (Log10)")+
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
        legend.position = "bottom") + facet_wrap(~Group)


pdf("HAP_lefse_output.pdf",width=16,height=11);
HAP_discriminant_contigs
dev.off()

# FIGURE 4E Chordiagram HAP

#Run MAASLIN2 on  counts of HAP-associated contigs with relative abundance of the core respiratory bacteriome

library(Maaslin2)

df_input_data = read.table(file = "counts_HAP_associated.txt",
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
col_fun = colorRamp2(c(-100,-50, -1 ,0, 1, 50), c("red", "pink", "white", "lightblue","blue"))

pdf("chordDiagram_HAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()

# FIGURE 4F Chordiagram NO HAP

data <- read.delim("Correlations_viral_Core_Bacteriome_NOHAP.txt", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(circlize)
library(scales)
library(reshape)

grid.col = c(Fusobacterium = "red", Haemophilus = "orange", Prevotella = "green", Streptococcus = "blue", Veillonella = "purple")
col_fun = colorRamp2(c(0, 1, 50,100), c("white", "lightblue","blue"))

pdf("chordDiagram_noHAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()
