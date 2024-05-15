############################################################################################################################
#### CODE USED TO GENERATE FIGURE 2
############################################################################################################################

# Alphadiv metrics :

library(vegan)

data<-read.delim("counts.txt", row.names = 1)
dataTransposed<-t(data)
Shannon<-diversity(dataTransposed, index = 'shannon')
write.table(Shannon,"ShannonDiversity.txt", sep = '\t')

dataTransposed<-t(data)
richness <- specnumber(meio.data) 
write.table(richness,"Richness.txt", sep = '\t')

# Alphadiv FIGURE 2A - Alphadiv HAP/noHAP

#Dot plot - Shannon

data<-read.delim("ShannonDiversity_over_time.txt", row.names = 1)
data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Shannon, type = "median_iqr")
stat.test <- data %>% 
  wilcox_test(Shannon ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Shannon ~ GROUP)

summary_stats <- data %>%
  group_by(GROUP) %>%
  summarise(median = median(Shannon),
            p25 = quantile(Shannon, 0.25),
            p75 = quantile(Shannon, 0.75))


Shannon_dotplot <- ggplot(data, aes(GROUP, Shannon)) +
  geom_dotplot(method = "histodot", binaxis = "y", stackratio = 1,stackdir = "center", binpositions="all", binwidth = 0.2, dotsize = 1, aes(fill = GROUP, stroke = GROUP)) +
  scale_fill_manual(values = c("NO_HAP" = "blue", "HAP" = "red")) +
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 7, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Shannon") +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Shannon") +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20), legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = FALSE),
    y = "Shannon index"
  ) +
  labs(
    subtitle = get_test_label(stat.test, detailed = FALSE),
    y = "Shannon index"
  )
Shannon_dotplot

pdf("Shannon_dotplot.pdf",width=12,height=8);
Shannon_dotplot
dev.off()


# Alphadiv FIGURE 2B - Alphadiv HAP/noHAP


#Dot plot - Richness

data<-read.delim("Richness_over_time.txt", row.names = 1)
data %>% sample_n_by(GROUP, size = 2)
data %>%
  group_by(GROUP) %>%
  get_summary_stats(Richness, type = "median_iqr")
stat.test <- data %>% 
  wilcox_test(Richness ~ GROUP) %>%
  add_significance()
stat.test
data %>% wilcox_effsize(Richness ~ GROUP)

summary_stats <- data %>%
  group_by(GROUP) %>%
  summarise(median = median(Richness),
            p25 = quantile(Richness, 0.25),
            p75 = quantile(Richness, 0.75))


Richness_dotplot <- ggplot(data, aes(GROUP, Richness)) +
  geom_dotplot(method = "histodot", binaxis = "y", stackratio = 1,stackdir = "center", binpositions="all", binwidth = 80, dotsize = 1, aes(fill = GROUP, stroke = GROUP)) +
  scale_fill_manual(values = c("NO_HAP" = "blue", "HAP" = "red")) + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 3000, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Richness") +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Richness") +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size = 20), legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = FALSE),
    y = "Richness"
  ) +
  labs(
    subtitle = get_test_label(stat.test, detailed = FALSE),
    y = "Richness"
  )
Richness_dotplot

pdf("Richness_dotplot.pdf",width=12,height=8);
Richness_dotplot
dev.off()


# PCOA FIGURE 2C - PCOA HAP/noHAP

(We followed the PCOA scripts from Montassier et al., Nat Med. 2023;29(11):2793-2804.)

data1<-read.delim("counts.txt", row.names = 1)
x<-t(data1)
m <- read.delim("metadata.txt",row.names = 1)
library(ggplot2)
library(vegan)
library(ape)
library(cowplot)
beta_table <- as.matrix(vegdist(x), method = "bray", na.rm = F)
PCOA <- pcoa(beta_table)$vectors
var_exp <- pcoa(beta_table)$values
# Run stats for differentiation centroids
beta_dist = as.dist(beta_table)
length(beta_dist)
# Run PERMANOVA
ad = adonis(beta_dist ~ m$HAP_condition, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
# Run stats for differentiation dispersion
beta_out <- betadisper(beta_dist, m$HAP_condition)
p_val_disp <- permutest(beta_out)$tab[1, 6]

# Run stats on the coordinates
PCOA < - PCOA[rownames(m),]
wilcox.test(PCOA[,1] ~ m$HAP_condition)
wilcox.test(PCOA[,2] ~ m$HAP_condition)
# Plot beta diversity PCoA
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep = "")
}
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
m <- cbind(m, rownames(PCOA))
m <- data.frame(lapply(m, as.character), stringsAsFactors=FALSE)
colnames(m)[ncol(m)] <- "SampleID"
PCOA <- merge(PCOA, m, by = "SampleID")
PCOA$PC1 <- as.numeric(as.character(PCOA$PC1))
PCOA$PC2 <- as.numeric(as.character(PCOA$PC2))
PCOA$PC3 <- as.numeric(as.character(PCOA$PC3))
PCOA$PC4 <- as.numeric(as.character(PCOA$PC4))
# Make PCoA plot
body_cols=c("HAP" = "red", "NO_HAP" = "blue")
body_PCOA <- ggplot(PCOA) +
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2",
                                              color = "HAP_condition")) +
  scale_color_manual(values=body_cols) +
  theme_cowplot(font_size = 7) +
  guides(color=F) +
  annotate("text", x = -0.45, y = -0.2, label= paste("P = ", p_val), size=2) +
  annotate("text", x = -0.45, y = -0.25, label= paste("R2 = ", round(r_sq,
                                                                     digits=3)), size=2) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(color=NA), axis.text.y =
          element_text(color=NA))

body_PCOA <- body_PCOA +
  stat_ellipse(aes(x = PC1, y = PC2, color = HAP_condition))

# Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = factor(PCOA$HAP_condition, levels=c("HAP",
                                                                   "NO_HAP")), y = "PC1", fill = "HAP_condition")) +
  scale_fill_manual(values=body_cols) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x = "", y= paste("PC1 (", round(var_exp$Relative_eig[1],digits=3)*100, "%)", sep = ""))
PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x =factor(PCOA$HAP_condition, levels=c("HAP",
                                                                  "NO_HAP")), y = "PC2", fill = "HAP_condition")) +
  scale_fill_manual(values=body_cols) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x ="", y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=
                                        3)*100, "%)", sep = "")) +
  theme(axis.text.x = element_text(color=NA))
# Compile the PCoA and boxes
top2 <- plot_grid(PC2_boxes, body_PCOA, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.3, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))
pdf("PCOA_HAP_noHAP.pdf",width=7,height=3.5);
together2
dev.off()

# NMDS FIGURE 2D - NMDS HAP/noHAP

library(vegan)
library(tidyverse)
library(ggplot2)

set.seed(1)
data1<-read.delim("counts.txt", row.names = 1)
dataTransposed1<-t(data1)
dist.1 <- vegdist(dataTransposed1, method = "bray")
metadata <- read.delim("metadata.txt")

ano = anosim(dataTransposed1, metadata$HAP_condition, distance = "bray", permutations = 9999)
ano
plot(ano)

nmds = metaMDS(dataTransposed1, distance = "bray")
nmds
plot(nmds)


ta.scores = as.data.frame(scores(nmds)$sites)
metadata$HAP_condition = metadata$HAP_condition

NMDS_HAP_noHAP = ggplot(metadata, aes(x = ta.scores$NMDS1, y = ta.scores$NMDS2)) + 
  geom_point(aes(size = 10, colour = GROUP), alpha = 0.3)+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "GROUP", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = c("RED", "BLUE","pink")) + stat_ellipse(geom = "polygon", aes(group = GROUP, color = GROUP, fill = GROUP), alpha = 0.05) +
  annotate("text", x = -2, y = 7, label = paste0("Stress: ", format(nmds$stress, digits = 4)), hjust = 2) +
  annotate("text", x = -2, y = 6, label = paste0("P=", format(ano$signif, digits = 4)), hjust = 3) +
  scale_x_continuous(breaks = c(-6,-3,0, 3, 6)) +
  scale_y_continuous(breaks = c(-6,-3,0, 3, 6))

NMDS_HAP_noHAP

pdf("NMDS_HAP_noHAP.pdf",width=10,height=5);
NMDS_HAP_noHAP
dev.off()


# FIGURE 2E - WBC

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_between_HAP_NOHAP.txt", stringsAsFactors = FALSE)
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

WBC_Between_violin_HAP_noHAP <- ggplot(data, aes(GROUP, WBC)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Between-NO HAP patients" = "blue"), name = "GROUP") + 
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

pdf("WBC_Between_violin_HAP_noHAP.pdf",width=12,height=9);
WBC_Between_violin_HAP_noHAP
dev.off()

# FIGURE 2E - Hellinger

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_between_HAP_NOHAP.txt", stringsAsFactors = FALSE)
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

Hellinger_Between_violin_HAP_noHAP <- ggplot(data, aes(GROUP, Hellinger)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Between-NO HAP patients" = "blue"), name = "GROUP") + 
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

pdf("Hellinger_Between_violin_HAP_noHAP.pdf",width=12,height=9);
Hellinger_Between_violin_HAP_noHAP
dev.off()

# FIGURE 2E - Sorensen

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_between_HAP_NOHAP.txt", stringsAsFactors = FALSE)
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

Sorensen_Between_violin_HAP_noHAP <- ggplot(data, aes(GROUP, Sorensen)) +
  geom_violin(aes(fill = GROUP), color = "black", trim = TRUE) +  
  geom_boxplot(aes(fill = GROUP), width = 0.05, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Between-HAP patients" = "red", "Between-NO HAP patients" = "blue"), name = "GROUP") + 
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

pdf("Sorensen_Between_violin_HAP_noHAP.pdf",width=12,height=9);
Sorensen_Between_violin_HAP_noHAP
dev.off()