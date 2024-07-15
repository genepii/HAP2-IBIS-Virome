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

# rarefaction

library("phyloseq")

library("ggplot2")

otu = read.delim(file="dRPKM_final.txt",header=T,row.names=1,check.names=F) 
tax = read.delim(file="tax.txt",header=T,row.names=1,check.names=F) 
meta = read.delim(file="IBIS_clinical_final.txt",header=T,row.names=1,check.names=F) 
otu <- round(otu, digits = 0)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

physeq = phyloseq(OTU, TAX)

physeq = merge_phyloseq(physeq, sample_data(meta))

#### 

require(parallel)
options(mc.cores= 16)
require(vegan)
## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#####


p <- ggrare(physeq, color = "HAP_condition", se = FALSE) 

p <- p + 
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~Patient, scales = "free") + 
  theme_bw()

plot(p + facet_wrap(~N_Atlanrea,scales = "free") + theme_bw())
ggsave(file="Rarefaction_curve_sample_virome_ibis_red.pdf",width=50,height=30, limitsize= FALSE)



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

#Loess plot - Shannon index

data<-read.delim("ShannonDiversity_over_time_before_HAP_onset.txt", row.names = 1)
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
data<-read.delim("Richness_over_time_before_HAP_onset.txt", row.names = 1)
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

data <- read.delim("Unifrac_sliding_window.txt", stringsAsFactors = FALSE)

custom_colors <- c("upcoming HAP" = "pink", "no HAP" = "blue")

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


WUF_dynamics_stat <- WUF_dynamics + scale_y_continuous(breaks = c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("5-4", "4-3", "3-2", "2-1"))+
  geom_text(data = data, aes(x = DAY, y = 0.6, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -15.5)

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
############################################################################################################################################################
