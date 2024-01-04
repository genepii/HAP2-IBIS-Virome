############################################################################################################################
#### CODE USED TO GENERATE FIGURE 6
############################################################################################################################

# FIGURE 6A

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
pdf("ALR_heatmap_prevhap.pdf",width=20,height=14);
ALR_heatmap
dev.off()


# RELAB heatmap FIGURE 6A

library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)

matrix <- read.delim("Relab_family.txt", header = TRUE, sep = "\t")

melted_data <- melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Family","Sample", "Value")


virus_colors <- c("#450659FF","#460B5EFF","#471063FF","#481668FF","#481A6CFF","#481E70FF","#0000FF","#482374FF",
                  "#472C7AFF","#472F7EFF","#463480FF","#453882FF","#433D84FF","#424086FF","#482778FF","#3F4889FF","#3D4D8AFF",
                  "#3C508BFF","#3A538BFF","#38588CFF","#375B8DFF","#355F8DFF","#33628DFF","#31668EFF","#30698EFF","#2E6D8EFF",
                  "#2D708EFF","#2C738EFF","#2A768EFF","#999999")


RELAB_heatmap_prevhap <- ggplot(melted_data, aes(fill = Family, y = Value, x = Sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = virus_colors) +
  labs(title = "Relative abundance of viral families",
       x = "Samples",
       y = "Relative abundance") +
  theme_classic() +
  geom_col(colour = "black", stat = "identity") +
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_blank())

pdf("RELAB_heatmap_prevhap.pdf",width=20,height=15);
RELAB_heatmap_prevhap
dev.off()

# FIGURE 6B - Betadiv - Bray_Curtis

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
  labs(x = "Days before and after HAP onset (Sliding window 3d)", y = "Bray Curtis Dissimilarity")


WBC_dynamics_stat <- WBC_dynamics + scale_y_continuous(breaks = c(0.9,0.95,1)) +scale_x_continuous(breaks = c(1, 2, 3), labels = c("5-3", "3-0", "0-4"))+
  geom_text(data = data, aes(x = DAY, y = 0.9, label = ifelse(Pvalue < 0.0001, "****", ifelse(Pvalue < 0.001, "***", ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "ns"))))), color = "black", size = 16, vjust = -17.5)

pdf("WBC_dynamics_prevhap.pdf",width=15,height=10);
WBC_dynamics_stat
dev.off()

# FIGURE 6C - Lefse HAP contigs

#Use LEfSe to identify discriminant contigs in the HAP signature 6 days before the HAP onset.
#Input table : ALR/Presence absence counts

#In shell, run :
#format_input.py Contig_ALR_HAP.txt HAP_lefse_input.in -c 1 -u 2 -o 1000000
#run_lefse.py HAP_lefse_input.in HAP_lefse_output.res


library(ggplot2)

data <- read.delim("Lefse_prevhap_output.res")

HAP_signature_contigs <- ggplot(data, aes(reorder(Taxon, LDA), LDA, fill = Group)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("Contigs associated with 5-3 period days before HAP onset" = "#6666ff","Contigs associated with 3-0 period days before HAP onset" ="#999999")) +
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


pdf("HAP_signature_contigs_prevhap.pdf",width=16,height=11);
HAP_signature_contigs
dev.off()

# FIGURE 6D - Core respiratory bacteriome relative abundance 

# (16S data from PREVHAP external cohort : unpublished data)
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

pdf("loessPlot_Core_bacteriome_relab_prevhap.pdf",width=15,height=10);
loessPlot
dev.off()

# FIGURE 6E Heatmapp of the 13 HAP-associated strictly to HAP in the external cohort

library(pheatmap)

matrix <- read.delim("RPKMcounts_13.txt", header = TRUE, sep = "\t", row.names = 1)
matrix <- na.omit(matrix)
my_sample_col <- read.table("metadata.txt", header = TRUE, sep = "\t",row.names = 1)
transformed_matrix <- log10(matrix+1)

ann_colors = list(
  GROUP = c(HAP = "red", NO_HAP = "blue")
)

Heatmap_13_contigs <- pheatmap(
  transformed_matrix,
  color=colorRampPalette(c("white","lawn green", "yellow green","lime green","dark green","sea green","cyan","turquoise","dark turquoise","cadet blue","deep sky blue","dodger blue","medium blue","blue","dark blue","navy","purple"))(1000),
  cluster_cols = TRUE,
  clustering_distance_cols =  "binary",
  cluster_rows = TRUE,
  annotation_col = my_sample_col,
  legend = TRUE,
  fontsize_row = 20,  # Add column annotations
  annotation_legend = TRUE, show_colnames = F,labels_col=" ", fontsize_col = 10, angle_col = 90,
  annotation_colors = ann_colors,
  clustering_method = "ward.D2"
)

pdf("Heatmap_13_contigs_prevhap.pdf",width=20,height=15);
Heatmap_13_contigs
dev.off()

# FIGURE 6F ROC curve based on the FFT model of an example contig among the 13

library(FFTrees)
library(pROC)
library(caret)

ibis <- read.delim("Presence_absence_IBIS.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#The dataset to predict contains the variable to predict (here Has.HAP) as "NA"

prevhap_blind <- read.delim("Presence_absence_PREVHAP.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

tree <- FFTrees(formula = Has.HAP ~ ., 
                
                data = ibis,
                
                main = "HAP prediction", 
                
                decision.labels = c("No HAP", "HAP"))

hap_predicted = predict(
  object = tree,
  newdata = prevhap_blind,
  tree = 1,
  type = "class")

hap_predicted
$expected_value
$predicted_value

Confusion_Matrix <- confusionMatrix(data=$predicted_value, reference = $expected_value)

Confusion_Matrix

table($expected_value,$predicted_value)

# Create a ROC curve
roc_curve <- roc($expected_value, as.numeric($predicted_value), percent = TRUE)

pdf("ROC_curve_prevhap.pdf")
plot.roc(roc_curve,
         print.auc = TRUE,
         auc.polygon = TRUE,
         grid = c(0.1, 0.2),
         grid.col = c("green", "red"),
         max.auc.polygon = TRUE,
         auc.polygon.col = "lightblue",
         print.thres = TRUE,
         add = FALSE,
         main = "ROC Curve")
dev.off()

##############################################################