############################################################################################################################
#### CODE USED TO GENERATE Supplementary figures
############################################################################################################################

# Supplementary FIGURE 1

library(vegan)
set.seed(1)
data1<-read.delim("counts.txt", row.names = 1)
dataTransposed1<-t(data1)
dist.1 <- vegdist(dataTransposed1, method = "bray")
metadata <- read.delim("metadata.txt")

adonis.test1.cont.subjID <- adonis2(dist.1 ~ Samples, data = metadata, permutations = 999, strata = metadata$BATCH)
adonis.test1.cont.subjID

dispersion <- betadisper(dist.1, group=metadata$BATCH, type=c("median","centroid"))
permutest(dispersion)


colors <- c("1" = "red", "2" = "blue", "3" = "green", "4" = "gray")


pdf("batch_pcoa.pdf",width=7,height=7);

batch_pcoa <- plot(dispersion, hull=FALSE, ellipse=TRUE, label.cex = 3, axes = c(1,2),
     lty = "solid", lwd =3, segments = TRUE, seg.col = "black", label = FALSE, main = "Bray dissimilarity between sequencing batches",
     col = colors, cex = 3, sub = "")
legend("topright", legend = names(colors), fill = colors, title = "Batch")

dev.off()

pdf("batch_boxplot.pdf",width=7,height=7);
batch_boxplot <- boxplot(dispersion, ylab = "Distance of centroids", col=colors)

dev.off()


# Supplementary FIGURE 2

library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)

matrix <- read.delim("Relative_abundance.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Family","Sample", "Value")



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

# Supplementary FIGURE 3

#Loess plot - Shannon index

data<-read.delim("ShannonDiversity_over_time.txt", row.names = 1)
custom_colors <- c("HAP" = "#CC0033", "NO_HAP" = "#0000FF")


Shannon_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Shannon, color = HAP_condition, group = HAP_condition)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,2,4,6)
  )+
  stat_smooth(metho1 = "loess", formula = y ~ x, aes(fill = HAP_condition), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14)) +
  scale_color_manual(values = custom_colors) + 
  scale_fill_manual(values = custom_colors) +  
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),    
    axis.title = element_text(size = 20),       
    legend.text = element_text(size = 20),   
    legend.title = element_blank()     
  ) +
  xlab("Days before and after HAP onset") +  
  ylab("Shannon index")   


Shannon_loessPlot

pdf("Shannon_loessPlot.pdf",width=12,height=8);
Shannon_loessPlot
dev.off()



# Supplementary FIGURE 4

#Loess plot - Richness

data<-read.delim("Richness_over_time.txt", row.names = 1)
custom_colors <- c("HAP" = "#CC0033", "NO_HAP" = "#0000FF")

Richness_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Richness, color = HAP_condition, group = HAP_condition)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,500,1000,1500,2000)
  )+
  stat_smooth(metho1 = "loess", formula = y ~ x, aes(fill = HAP_condition), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14)) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20),      
    legend.text = element_text(size = 20),    
    legend.title = element_blank()     
  ) +
  xlab("Days before and after HAP onset") + 
  ylab("Richness")


Richness_loessPlot

pdf("Richness_loessPlot.pdf",width=12,height=8);
Richness_loessPlot
dev.off()


# Supplementary FIGURE 5

#Loess plot - Shannon index

data<-read.delim("IBIS_clinical_final_Upcoming_noHAP.txt", row.names = 1)
custom_colors <- c("Upcoming_HAP" = "pink", "NO_HAP" = "#0000FF")

# loess with custom colors
Shannon_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Shannon, color = GROUP, group = GROUP)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5))+
  stat_smooth(metho1 = "loess", formula = y ~ x, aes(fill = GROUP), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6,-5, -4, -3,-2,-1)) +
  scale_color_manual(values = custom_colors) +  # Add custom colors for points
  scale_fill_manual(values = custom_colors) +   # Add custom colors for stat_smooth
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),        # Adjust size of axis text
    axis.title = element_text(size = 20),       # Adjust size of axis titles
    legend.text = element_text(size = 20),      # Adjust size of legend text
    legend.title = element_blank()      # Adjust size of legend title
  ) +
  xlab("Days relative to HAP onset") +  # Add or modify the x-axis label
  ylab("Shannon index")    # Add or modify the y-axis label


Shannon_loessPlot

pdf("Shannon_loessPlot_up.pdf",width=12,height=8);
Shannon_loessPlot
dev.off()

# Supplementary FIGURE 6

#Loess plot - Richness
data<-read.delim("IBIS_clinical_final_Upcoming_noHAP.txt", row.names = 1)
custom_colors <- c("Upcoming_HAP" = "pink", "NO_HAP" = "#0000FF")

# loess with custom colors
Richness_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Richness, color = GROUP, group = GROUP)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000)
  )+
  stat_smooth(metho1 = "lm", formula = y ~ x, aes(fill = GROUP), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6,-5, -4, -3,-2,-1)) +
  scale_color_manual(values = custom_colors) +  # Add custom colors for points
  scale_fill_manual(values = custom_colors) +   # Add custom colors for stat_smooth
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),        # Adjust size of axis text
    axis.title = element_text(size = 20),       # Adjust size of axis titles
    legend.text = element_text(size = 20),      # Adjust size of legend text
    legend.title = element_blank()      # Adjust size of legend title
  ) +
  xlab("Days relative to HAP onset") +  # Add or modify the x-axis label
  ylab("Richness")    # Add or modify the y-axis label


Richness_loessPlot

pdf("Richness_loessPlot_up.pdf",width=12,height=8);
Richness_loessPlot
dev.off()



# Supplementary FIGURE 7

Dynamics Weighted Unifrac Distance (WUF) in 16S rRNA (We used 16S rRNA data from Montassier et al., Nat Med. 2023;29(11):2793-2804)

#Calculate Unifrac distance :

library(phyloseq)
library(ggplot2)
library(ape)

otu = read.delim(file="16S_count.txt",header=T,row.names=1,check.names=F) 
tax = read.delim(file="taxonomy.txt",header=T,row.names=1,check.names=F)
metadata = read.delim(file="metadata.txt",header=T,row.names=1,check.names=F)

colnames(tax) = c(
  "Superkingdom",
  "Phylum",
  "Class"	,
  "Order"	,
  "Family",
  "Genus"	,
  "Species",
  "Strain")

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))
physeq = phyloseq(OTU, TAX)
physeq = merge_phyloseq(physeq, sample_data(meta))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
dis<- UniFrac(physeq1, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
write.table(dis,"WeightedUnifrac.txt", sep = '\t')

#Plot Unifrac dynamics :

data <- read.delim("WUF_dynamics.txt", stringsAsFactors = FALSE)

custom_colors <- c("HAP" = "red", "NO_HAP" = "blue")

WUF_dynamics <- ggplot(data, aes(x = DAY, y = WUF, group = GROUP, color = GROUP)) +
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
  labs(x = "Days before HAP onset (Sliding window 2d)", y = "Weighted Unifrac Distance")


WUF_dynamics_stat <- WUF_dynamics + scale_y_continuous(breaks = c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("5-3", "4-2", "3-1", "2-0"))+
  geom_text(data = data, aes(x = DAY, y = 0.5, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -17.3)

pdf("WUF_dynamics.pdf",width=15,height=10);
WUF_dynamics_stat
dev.off()

# Supplementary FIGURE 8


library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)

matrix <- read.delim("RELAB_prevhap_fam_agg.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Virus","Sample", "Value")

virus_colors <- c("#440154FF","#450559FF","#460A5DFF","#470F62FF","#481467FF","#48186AFF","#481C6EFF","#482072FF",
                  "#482576FF","#482979FF","#472D7BFF","#46307EFF","#463480FF","#453883FF","#0000FF","#424086FF",
                  "#414487FF","#3F4788FF","#3E4B89FF","#3C4F8AFF","#3B528BFF","#39558CFF","#38598CFF","#365C8DFF",
                  "#35608DFF","#33638DFF","#31668EFF","#30698EFF","#2F6C8EFF","#2E6F8EFF","#2C728EFF","#2B758EFF",
                  "#2A788EFF","#297B8EFF","#277E8EFF","#26828EFF","#25848EFF","#24878EFF","#238A8DFF","#228D8DFF",
                  "#21908CFF","#20938CFF","#1F968BFF","#1F998AFF","#1E9C89FF","#1F9F88FF","#1FA287FF","#21A585FF",
                  "#22A884FF","#25AB82FF","#27AD81FF","#2BB17FFF","#2FB47CFF","#34B679FF","#38B977FF","#3EBC74FF",
                  "#43BF71FF","#49C16DFF","#50C46AFF","#56C667FF","#5DC863FF","#64CB5FFF","#6BCD5AFF","#72D056FF",
                  "#7AD151FF","#82D34DFF","#8AD547FF","#92D742FF","#9AD93CFF","#A2DA37FF","#AADC32FF","#B3DD2CFF",
                  "#BBDF27FF","#C4E022FF","#CDE11DFF","#999999","#DEE318FF","#E6E419FF","#EEE51CFF","#F6E620FF",
                  "#FDE725FF")



RELAB_heatmap_prevhap <- ggplot(melted_data, aes(fill = Virus, y = Value, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = virus_colors) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  theme(legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_blank())

pdf("RELAB_heatmap_prevhap.pdf",width=20,height=15);
RELAB_heatmap_prevhap
dev.off()

# Supplementary FIGURE 9

#wbc boxplot

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_53.txt", stringsAsFactors = FALSE)
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

WBC_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), WBC)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Bray-Curtis dissimilarity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none")


pdf("WBC_boxplot.pdf",width=12,height=9);
WBC_boxplot
dev.off()

# hellinger boxplot

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_53.txt", stringsAsFactors = FALSE)
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

Hellinger_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), Hellinger)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 1.4, bracket.size = 2) +
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

pdf("Hellinger_boxplot.pdf",width=12,height=9);
Hellinger_boxplot
dev.off()

#sorensen boxplot 

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_53.txt", stringsAsFactors = FALSE)
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

Sorensen_boxplot <- ggplot(data, aes(reorder(GROUP, GROUP, function(x) -sum(x == "Upcoming HAP")), Sorensen)) +
  geom_boxplot(aes(fill = GROUP), width = 2, color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c("Upcoming HAP" = "pink", "NO HAP" = "blue"), name = "GROUP") +  
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, y.position = 1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE), y = "Sorensen dissimilarity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 30),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "none")


pdf("Sorensen_boxplot.pdf",width=12,height=9);
Sorensen_boxplot
dev.off()

############################################################################################################################################################
