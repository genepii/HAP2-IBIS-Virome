############################################################################################################################
#### CODE USED TO GENERATE Supplementary figures
############################################################################################################################

# Supplementary FIGURE 1

library(vegan)
set.seed(1)
data1<-read.delim("RPKMcounts.txt", row.names = 1)
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
library(reshape)
library(magrittr)

matrix <- read.delim("dRPKM_final_family.txt", header = TRUE, sep = "\t", row.names = 1)
Relab <- matrix %>%
  apply(2, function(x) x/sum(x)) %>%
  `*`(100)

head(Relab)

write.table(Relab,file="Relab_family.txt",quote=F,sep="\t",col.names=NA)

matrix <- read.delim("Relab_family.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Family","Sample", "Value")



virus_colors <- c("#450659FF","#460B5EFF","#471063FF","#481668FF","#481A6CFF","#481E70FF","#482374FF","#482778FF",
                  "#472C7AFF","#472F7EFF","#463480FF","#453882FF","#433D84FF","#0000FF","#424086FF","#3F4889FF","#3D4D8AFF",
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
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
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

# PCOA HAP/noHAP

(We followed the PCOA scripts from Montassier et al., Nat Med. 2023;29(11):2793-2804.)
set.seed(1)
data1<-read.delim("RPKMcounts.txt", row.names = 1)
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
ad = adonis(beta_dist ~ m$GROUP, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
# Run stats for differentiation dispersion
beta_out <- betadisper(beta_dist, m$GROUP)
p_val_disp <- permutest(beta_out)$tab[1, 6]

# Run stats on the coordinates
PCOA < - PCOA[rownames(m),]
wilcox.test(PCOA[,1] ~ m$GROUP)
wilcox.test(PCOA[,2] ~ m$GROUP)
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
body_cols=c("HAP" = "red", "no_HAP" = "blue")
body_PCOA <- ggplot(PCOA) +
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2",
                                              color = "GROUP")) +
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
  stat_ellipse(aes(x = PC1, y = PC2, color = GROUP))

# Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = factor(PCOA$GROUP, levels=c("HAP",
                                                                   "no_HAP")), y = "PC1", fill = "GROUP")) +
  scale_fill_manual(values=body_cols) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x = "", y= paste("PC1 (", round(var_exp$Relative_eig[1],digits=3)*100, "%)", sep = ""))
PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x =factor(PCOA$GROUP, levels=c("HAP",
                                                                  "no_HAP")), y = "PC2", fill = "GROUP")) +
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

# Supplementary FIGURE 6

# NMDS HAP/noHAP

library(vegan)
library(tidyverse)
library(ggplot2)

set.seed(1)
data1<-read.delim("RPKMcounts.txt", row.names = 1)
dataTransposed1<-t(data1)
dist.1 <- vegdist(dataTransposed1, method = "bray")
metadata <- read.delim("metadata.txt")

ano = anosim(dataTransposed1, metadata$GROUP, distance = "bray", permutations = 9999)
ano
plot(ano)

nmds = metaMDS(dataTransposed1, distance = "bray")
nmds
plot(nmds)


ta.scores = as.data.frame(scores(nmds)$sites)
metadata$HAP_condition = metadata$GROUP

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
  annotate("text", x = -2, y = 6, label = paste0("P=", format(ano$signif, digits = 4)), hjust = 3)+
  scale_x_continuous(breaks = c(-8,-4,0, 4)) +
  scale_y_continuous(breaks = c(-8,-4,0, 4))

NMDS_HAP_noHAP

pdf("NMDS_HAP_noHAP.pdf",width=10,height=8);
NMDS_HAP_noHAP
dev.off()

# Supplementary FIGURE 7

#Loess plot - Shannon index

data<-read.delim("ShannonDiversity_over_time_before_HAP_onset.txt", row.names = 1)
custom_colors <- c("Upcoming_HAP" = "pink", "NO_HAP" = "#0000FF")

# loess with custom colors
Shannon_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Shannon, color = GROUP, group = GROUP)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5))+
  stat_smooth(metho1 = "loess", formula = y ~ x, aes(fill = GROUP), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6,-5, -4, -3,-2,-1)) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),      
    axis.title = element_text(size = 20),       
    legend.text = element_text(size = 20),   
    legend.title = element_blank()  
  ) +
  xlab("Days relative to HAP onset") + 
  ylab("Shannon index")    


Shannon_loessPlot

pdf("Shannon_loessPlot_up.pdf",width=12,height=8);
Shannon_loessPlot
dev.off()

# Supplementary FIGURE 8

#Loess plot - Richness
data<-read.delim("Richness_over_time_before_HAP_onset.txt", row.names = 1)
custom_colors <- c("Upcoming_HAP" = "pink", "NO_HAP" = "#0000FF")

# loess with custom colors
Richness_loessPlot <- ggplot(data, aes(x = HAP_onset, y = Richness, color = GROUP, group = GROUP)) +
  geom_point(shape =20, size=7) +
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000)
  )+
  stat_smooth(metho1 = "lm", formula = y ~ x, aes(fill = GROUP), alpha = 0.3) +
  scale_x_continuous(breaks = c(-6,-5, -4, -3,-2,-1)) +
  scale_color_manual(values = custom_colors) + 
  scale_fill_manual(values = custom_colors) +   
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),      
    axis.title = element_text(size = 20),      
    legend.text = element_text(size = 20),   
    legend.title = element_blank()     
  ) +
  xlab("Days relative to HAP onset") +  
  ylab("Richness")   


Richness_loessPlot

pdf("Richness_loessPlot_up.pdf",width=12,height=8);
Richness_loessPlot
dev.off()


# Supplementary FIGURE 9

library(tidyverse)

# define data from correlations test

data_long <- data %>%
  pivot_longer(cols = c(Negative, Positive), names_to = "Interaction", values_to = "Count")


data_long <- data_long %>%
  group_by(Group) %>%
  mutate(Total = sum(Count))

data_long <- data_long %>%
  mutate(Percentage = Count / Total * 100)

library(ggplot2)

p <- ggplot(data_long, aes(x = "", y = Percentage, fill = Interaction)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste0(round(Percentage, 0), "%")),
            position = position_stack(vjust = 0.5),
            size = 20) +  
  scale_fill_manual(values = c("pink", "lightblue")) +
  facet_wrap(~ Group)  

p



pdf("pie.pdf",width=20,height=10)

p

dev.off()

# Supplementary FIGURE 10

library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)

matrix <- read.delim("RELAB_prevhap.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Virus","Sample", "Value")

virus_colors <- c("#440154FF","#450559FF","#460A5DFF","#470F62FF","#481467FF","#48186AFF","#481C6EFF","#482072FF",
                  "#482576FF","#482979FF","#472D7BFF","#46307EFF","#463480FF","#453883FF","#424086FF","#0000FF",
                  "#414487FF","#3F4788FF","#3E4B89FF","#3C4F8AFF","#3B528BFF","#39558CFF","#38598CFF","#365C8DFF",
                  "#35608DFF","#33638DFF","#31668EFF","#30698EFF","#2F6C8EFF","#2E6F8EFF","#2C728EFF","#2B758EFF",
                  "#2A788EFF","#297B8EFF","#277E8EFF","#26828EFF","#25848EFF","#24878EFF","#238A8DFF","#228D8DFF",
                  "#21908CFF","#20938CFF","#1F968BFF","#1F998AFF","#1E9C89FF","#1F9F88FF","#1FA287FF","#21A585FF",
                  "#22A884FF","#25AB82FF","#27AD81FF","#2BB17FFF","#2FB47CFF","#34B679FF","#38B977FF","#3EBC74FF",
                  "#43BF71FF","#49C16DFF","#50C46AFF","#56C667FF","#5DC863FF","#64CB5FFF","#6BCD5AFF","#72D056FF",
                  "#7AD151FF","#82D34DFF","#8AD547FF","#92D742FF","#9AD93CFF","#A2DA37FF","#AADC32FF","#B3DD2CFF",
                  "#BBDF27FF","#C4E022FF","#CDE11DFF","#00FF00FF","#DEE318FF","#E6E419FF","#EEE51CFF","#F6E620FF",
                  "#FDE725FF","#333333", "#999999", "#0000FFFF", "#FFFF00FF", "#00FFFFFF", "#FF00FFFF", "#FF7F00FF")

RELAB_heatmap_prevhap <- ggplot(melted_data, aes(fill = Virus, y = Value, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = virus_colors) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  theme(legend.text = element_text(size = 16),
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


# Supplementary FIGURE 11

Phage lifestyle Barplot

#### Vibrant results

#VIBRANT to predict lifestyles (lysogenic/lytic) for vOTUs with minimum sequence length of 1000bp and containing at least 4 ORFs (open readings frames) 

singularity shell vibrant.sif
VIBRANT_run.py -i viral_vOTUs.fasta -t 114 -folder VIBRANT_results -virome

library(ggplot2)

# Define the data (Vibrant output)

# Define custom virus colors
style_colors <- c("Lytic" = "wheat1", "Lysogenic" = "orange3")

# Create the plot
plot <- ggplot(data, aes(fill = Virus, y = Percent, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = style_colors) +
  labs(x = "Cohort",
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


pdf("Vibrant.pdf",width=15,height=20);
plot
dev.off()

# Supplementary FIGURE 12
library(tidyverse)

# define data from correlations test

data_long <- data %>%
  pivot_longer(cols = c(Negative, Positive), names_to = "Interaction", values_to = "Count")


data_long <- data_long %>%
  group_by(Group) %>%
  mutate(Total = sum(Count))

data_long <- data_long %>%
  mutate(Percentage = Count / Total * 100)

library(ggplot2)

p <- ggplot(data_long, aes(x = "", y = Percentage, fill = Interaction)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste0(round(Percentage, 0), "%")),
            position = position_stack(vjust = 0.5),
            size = 20) +  
  scale_fill_manual(values = c("pink", "lightblue")) +
  facet_wrap(~ Group)  

p



pdf("pie.pdf",width=20,height=10)

p

dev.off()

############################################################################################################################################################
