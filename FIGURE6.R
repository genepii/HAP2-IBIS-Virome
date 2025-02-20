############################################################################################################################
#### CODE USED TO GENERATE FIGURE 6
############################################################################################################################

# A

#Blast viral signature in validation cohort in shell
#blastn -num_threads 32 -query viral_signature.fasta -subject Virome_validation.fasta -task blastn -max_target_seqs 1 -evalue 1e-10 -max_hsps 1 -qcov_hsp_perc 50 -outfmt '6 std stitle qseq sseq' >> "blast.tsv"

# plot blast
library(ggplot2)
library(dplyr)
library(reshape2)

data <- read.table("blast.txt", header = TRUE, sep = "\t")

# Reshape data if necessary
data_melted <- data %>%
  mutate(ibis = as.factor(ibis)) %>%
  melt(id.vars = c("ibis", "prevhap"))

blast<-ggplot(data_melted, aes(x = reorder(prevhap, value), y = reorder(ibis, value), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_light() +
  labs(title = "Heatmap of Alignment and Identity",
       x = "PREVHAP",
       y = "IBIS",
       fill = "Value") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16))

pdf("blast.pdf", height = 5, width = 7)
blast
dev.off()

# B plot conserved signature features

library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)

data <- read.table("Conserved_Signature.txt", header = TRUE, sep = "\t")

plot_sequence_length <- ggplot(data, aes(x = reorder(Origin, Origin_len), y = Origin_length)) +
  geom_bar(stat = "identity", fill = "purple4") +
  coord_flip() +
  labs(
    title = "vOTU Lengths",
    x = "vOTUs",
    y = "Length (bp)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    legend.text = element_text(size = 16)
  )

plot_protein_length <- ggplot(data, aes(x = reorder(Sequence_Name, Sequence_Length), y = Sequence_Length)) +
  geom_bar(stat = "identity", fill = "gold") +
  coord_flip() +
  labs(
    title = "Protein Lengths",
    x = "Identified proteins",
    y = "Length (bp)"
  ) +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    legend.text = element_text(size = 16)
  )

lifestyle_data <- data %>%
  count(Lifestyle) %>%
  mutate(Percentage = n / sum(n) * 100) 

colnames(lifestyle_data) <- c("Lifestyle", "Count", "Percentage")

plot_lifestyle <- ggplot(lifestyle_data, aes(x = Lifestyle, y = Percentage, fill = Lifestyle)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "D") +
  labs(
    title = "",
    x = "Lifestyle",
    y = "Percentage (%)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

tax_data <- data.frame(table(data$Taxonomical_Information))
colnames(tax_data) <- c("Taxon", "Count")

plot_taxonomy <- ggplot(tax_data, aes(x = "", y = Count, fill = Taxon)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_brewer(palette = "Set3") +
  coord_polar("y") +
  labs(
    title = "Taxonomical composition",
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 16))

prot_data <- data.frame(table(data$Protein_Information))
colnames(prot_data) <- c("Protein", "Count")

plot_proteins <- ggplot(prot_data, aes(x = "", y = Protein, fill = Protein)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_brewer(palette = "Set3") +
  coord_polar("y") +
  labs(
    title = "Protein composition",
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 16))

host_data <- data.frame(table(data$Host))
colnames(host_data) <- c("Host", "Count")

plot_host <- ggplot(host_data, aes(x = "", y = Count, fill = Host)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_viridis_d(option = "C") +
  coord_polar("y") +
  labs(
    title = "Host composition",
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 16))
		
top <- plot_grid(plot_sequence_length, plot_taxonomy, plot_protein_length,
                  plot_host,plot_lifestyle, plot_proteins, nrow =3, ncol=2, rel_widths =c(0.3,0.7))


pdf("plot.pdf", height = 20, width = 20)
top
dev.off()

# C PLSDA

library(mixOmics)


Clinical <- read.delim("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
Y <- factor(Clinical$GROUP)
Y_numeric <- as.numeric(factor(Y))
summary(Y)
X1 <- read.delim("dRPKM.txt", header = TRUE, sep = "\t", row.names = 1)
X1<-t(X1)

legend.color <- ifelse(Clinical$GROUP == "HAP", "blue", "PINK")
result.plsda.srbct <- plsda(X1, Y)
plsda.res <- plsda(X1, Y, ncomp = 10)
background.max <- background.predict(plsda.res, 
                                     comp.predicted = 2,
                                     dist = 'max.dist')

pdf("plsda.pdf", width = 5, height = 5)
plotIndiv(plsda.res, comp = 1:2, group = Clinical$GROUP,star = T, style = "lattice",
          col.per.group = c("HAP" ="pink", "NOHAP" = "blue"),
          ind.names = FALSE, title = 'Maximum distance',
          legend = TRUE,  background = background.max)
dev.off()

# D FFTREES And ten-fold cross-validations

library(pROC)
library(FFTrees)
library(caret) 

df <- read.delim("data.txt", header = TRUE, stringsAsFactors = FALSE,row.names = 1)
df$Has.HAP <- as.logical(df$Has.HAP)
folds <- createFolds(df$Has.HAP, k = 10, list = TRUE, returnTrain = FALSE)
all_pred_probs <- vector("list", length = 10)
# Run 10-fold cross-validation
for (i in 1:10) {
  test_indices <- folds[[i]]
  train_data <- df[-test_indices, ]
  test_data <- df[test_indices, ]

  train_data$Has.HAP <- as.logical(train_data$Has.HAP)
  test_data$Has.HAP <- as.logical(test_data$Has.HAP)
  
  custom.fft_model <- FFTrees(
    formula = Has.HAP ~ ., 
    data = train_data, 
    decision.labels = c("No HAP", "HAP") 
  )
  
  # Predict probabilities for the test set
  pred_prob <- predict(custom.fft_model, newdata = test_data, type = "prob")
  
  if (!is.null(pred_prob)) {
    all_pred_probs[[i]] <- pred_prob[, 2]
  } else {
    warning(paste("Prediction failed for fold", i))
  }
}

pred_probs_combined <- unlist(all_pred_probs)

if (length(pred_probs_combined) > 0) {
  true_labels <- df$Has.HAP
  roc_obj_cv <- roc(true_labels, pred_probs_combined, percent = TRUE)
  ci_auc <- ci.auc(roc_obj_cv)
  pdf("roc.pdf", height = 5, width = 5)
  plot.roc(roc_obj_cv,
           print.auc = TRUE,
           auc.polygon = TRUE,
           grid = c(0.1, 0.2),
           grid.col = c('green', 'red'),
           max.auc.polygon = TRUE,
           auc.polygon.col = 'white',
           print.thres = TRUE,
           col = "blue",
           main = 'ROC Curve with 10-Fold CV (FFTrees)')
  abline(a = 100, b = -1, col = "red", lty = 2)
  dev.off()
} else {
  warning("Predictions are empty, please check the model training and prediction steps.")
}

# E
#Run TKNA python script and plos using Cytoscape

# F plot TKNA bibc

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)
library(ggrepel)

data <- read.csv2("probabilities_to_randomly_find_nodes.csv")
data <- data %>%
  mutate(
    Probability = as.numeric((Probability)) * 100,
    Observed_BiBC = as.numeric(Observed_BiBC )
  )

p <- ggplot(data, aes(x = Observed_BiBC, y = Observed_degree, color = Probability)) +
  geom_point(size = 10, alpha = 0.3) +  # Points
  geom_text_repel(aes(label = Node), size = 5, max.overlaps = 10) + 
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +  # Minimal theme
  labs(
    title = "",
    x = "Observed BiBC",
    y = "Observed Degree",
    color = "Likelihood of random
  finding probability"
  ) +
  theme(legend.position = "right")

print(p)

pdf("bibc.pdf", height = 7, width = 7)
p
dev.off()

##############################################################