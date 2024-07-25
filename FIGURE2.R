############################################################################################################################
#### CODE USED TO GENERATE FIGURE 2
############################################################################################################################

# Alphadiv metrics :

library(vegan)

data<-read.delim("RPKMcounts.txt", row.names = 1)
dataTransposed<-t(data)
Shannon<-diversity(dataTransposed, index = 'shannon')
write.table(Shannon,"ShannonDiversity.txt", sep = '\t')

dataTransposed<-t(data)
richness <- specnumber(dataTransposed) 
write.table(richness,"Richness.txt", sep = '\t')

# Alphadiv FIGURE 2A - Alphadiv HAP/noHAP

#Dot plot - Shannon

data<-read.delim("diversity_metadata.txt", row.names = 1)
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
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40), legend.position = "none") +
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

data<-read.delim("diversity_metadata.txt", row.names = 1)
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
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40), legend.position = "none") +
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



# FIGURE 2C - WBC

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("WBC_between.txt", stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Inter HAP" = "red", "Inter no HAP" = "blue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Bray-Curtis dissimilarity"
  )

pdf("WBC_Between_violin_HAP_noHAP.pdf",width=12,height=9);
WBC_Between_violin_HAP_noHAP
dev.off()

# FIGURE 2C - Hellinger

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Hellinger_between.txt", stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Inter HAP" = "red", "Inter no HAP" = "blue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.5, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Hellinger dissimilarity"
  )

pdf("Hellinger_Between_violin_HAP_noHAP.pdf",width=12,height=9);
Hellinger_Between_violin_HAP_noHAP
dev.off()

# FIGURE 2C - Sorensen

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read.delim("Sorensen_between.txt", stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Inter HAP" = "red", "Inter no HAP" = "blue"), name = "GROUP") + 
  stat_pvalue_manual(stat.test, tip.length = 0, size = 20, y.position = 1.1, bracket.size = 2) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 40),
        plot.subtitle = element_text(size = 30),
        legend.text = element_text(size = 40),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(stat.test, detailed = TRUE),
    y = "Sorensen dissimilarity"
  )

pdf("Sorensen_Between_violin_HAP_noHAP.pdf",width=12,height=9);
Sorensen_Between_violin_HAP_noHAP
dev.off()
