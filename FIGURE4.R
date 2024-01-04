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

custom_colors <- c("HAP" = "red", "NO_HAP" = "blue")

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

custom_colors <- c("HAP" = "red", "NO_HAP" = "blue")

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

custom_colors <- c("HAP" = "red", "NO_HAP" = "blue")

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

custom_colors <- c("HAP" = "red", "NO_HAP" = "blue")

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

# FIGURE 4C - Lefse HAP contigs

#Use LEfSe to identify discriminant contigs in the HAP signature 6 days before the HAP onset.
#Input table : ALR/Presence absence counts

#In shell, run :
#format_input.py Contig_ALR_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("HAP_lefse_output.res")

HAP_signature_contigs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Contigs associated with 5-3 period days before HAP onset" = "#cc00cc","Contigs associated with 3-0 period days before HAP onset" ="#ffcc00")) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "LEfSe of contigs in HAP signature", x = "HAP patients discriminant viral contigs", y = "LDA score (Log10)")+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        plot.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "bottom") 


pdf("HAP_signature_contigs.pdf",width=16,height=11);
HAP_signature_contigs
dev.off()
