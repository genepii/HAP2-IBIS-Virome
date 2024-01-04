############################################################################################################################
#### CODE USED TO GENERATE FIGURE 5
############################################################################################################################

# FIGURE 5A Chordiagram HAP

#Run MAASLIN2 on RPKM counts of HAP-associated contigs/NO HAP associated contigs with relative abundance of the core respiratory bacteriome

library(Maaslin2)

df_input_data = read.table(file = "RPKMcounts_associated_HAP/NO HAP.txt",
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
col_fun = colorRamp2(c(-50, -1 ,0, 1, 50), c("red", "pink", "white", "lightblue","blue"))

pdf("chordDiagram_HAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()

# FIGURE 5B Chordiagram NO HAP

data <- read.delim("Correlations_viral_Core_Bacteriome_NOHAP.txt", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)
library(circlize)
library(scales)
library(reshape)

grid.col = c(Fusobacterium = "red", Haemophilus = "orange", Prevotella = "green", Streptococcus = "blue", Veillonella = "purple")
col_fun = colorRamp2(c(-50, -1 ,0, 1, 50), c("red", "pink", "white", "lightblue","blue"))

pdf("chordDiagram_noHAP.pdf",width=20,height=20)
chordDiagram(data, big.gap = 10, symmetric = TRUE, annotationTrack = c("grid","name"), col = col_fun, grid.col = grid.col, scale = FALSE)
dev.off()


# FIGURE 5C - Dynamics Weighted Unifrac Distance (WUF) in 16S (We used 16s data from Montassier et al., Nat Med. 2023;29(11):2793-2804)

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

# FIGURE 5D - Core respiratory bacteriome relative abundance 

library(ggplot2)

data <- read.delim("Core_bacteriome_relative_abundance.txt")


loessPlot <- ggplot(data, aes(x = time, y = Value, color = Bacteria)) +
  geom_point() +
  stat_smooth(method = "loess", formula = y ~ x, aes(fill = Bacteria), alpha = 0.3) + 
  scale_x_continuous(breaks = c(-6, -5, -4, -3, -2, -1, 0)) +  
  theme_classic() +
  facet_grid(~GROUP, scales = "free") +
  theme(
    axis.text = element_text(size = 20),        
    axis.title = element_text(size = 20),     
    legend.text = element_text(size = 20),     
    legend.title = element_blank(),     
    strip.text = element_text(size = 20)        
  ) +
  xlab("Days before HAP onset") + 
  ylab("Respiratory Core Microbiome Relative abundance (%)") +
  coord_cartesian(ylim = c(0, NA)) 

loessPlot

pdf("loessPlot_Core_bacteriome_relab_ibis.pdf",width=15,height=10);
loessPlot
dev.off()