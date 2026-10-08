cat("\nDone! Results saved to Desktop!\n")
mydata <- read.csv("C:/CDS FASTA/COMBINED RESULTS.csv")
mydata <- mydata[, names(mydata) != ""]
cat("=== ONE-WAY ANOVA: ENC by Virus Type ===\n\n")
anova_result <- aov(ENC ~ Type, data = mydata)
print(summary(anova_result))
cat("\nMean ENC per group:\n")
print(tapply(mydata$ENC, mydata$Type, function(x)
  c(Mean=round(mean(x),4), SD=round(sd(x),4),
    Min=round(min(x),4), Max=round(max(x),4))))
cat("\n=== PEARSON CORRELATIONS: ENC vs GC Content ===\n\n")
pairs_list <- list(
  "GC%"   = "GC",
  "GC1%"  = "GC1",
  "GC2%"  = "GC2",
  "GC3s%" = "GC3s",
  "GC12%" = "GC12"
)
results_df <- data.frame()
for (name in names(pairs_list)) {
  col <- pairs_list[[name]]
  test <- cor.test(mydata$ENC, mydata[[col]], method = "pearson")
  cat(sprintf("ENC vs %s: r = %.4f, p = %.6f %s\n",
    name, test$estimate, test$p.value,
    ifelse(test$p.value < 0.05, "(SIGNIFICANT)", "(not significant)")))
  results_df <- rbind(results_df, data.frame(
    Comparison  = paste("ENC vs", name),
    r_value     = round(test$estimate, 4),
    p_value     = round(test$p.value, 6),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No")))
}
write.csv(results_df, "C:/Users/hp/Documents/Pearson_correlations.csv", row.names=FALSE)
cat("\nDone! File saved to Documents folder!\n")
cat("Go to: C:/Users/hp/Documents/ to find Pearson_correlations.csv\n")
q()
# ============================================================
# ENC-Plot Analysis 
# ============================================================================
# GENE-LEVEL COMPARISON: hRSV-A vs hRSV-B
# ============================================================================
library(tidyverse)
# Load data 
— use forward slashes even on Windows
data <- read.csv("C:/Users/hp/Downloads/COMBINED RESULTS (2).csv", stringsAsFactors = FALSE)
# Extract gene name from title column (before first underscore)
data$Gene <- sapply(strsplit(data$title, "_"), function(x) x[1])
# Metrics you want to compare gene-wise
metrics <- c("ENC", "CAI", "GC3s", "GC", "Fop", "CBI")
genes <- unique(data$Gene)
# ============================================================================
# GENE-WISE ANOVA TABLE (mimics your screenshot format)
# ============================================================================
results <- data.frame()
for (gene in genes) {
  gene_data <- data %>% filter(Gene == gene)
  a_vals <- gene_data %>% filter(Type == "hRSV-A")
  b_vals <- gene_data %>% filter(Type == "hRSV-B")
  for (metric in metrics) {
    a <- a_vals[[metric]]
    b <- b_vals[[metric]]
    # One-way ANOVA (equivalent to t-test for 2 groups, gives F and p)
    model <- aov(gene_data[[metric]] ~ gene_data$Type)
    f_val <- summary(model)[[1]]$'F value'[1]
    p_val <- summary(model)[[1]]$'Pr(>F)'[1]
    results <- rbind(results, data.frame(
      Gene = gene,
      Metric = metric,
      hRSV_A_Mean = round(mean(a, na.rm = TRUE), 5),
      hRSV_A_SD = round(sd(a, na.rm = TRUE), 6),
      hRSV_B_Mean = round(mean(b, na.rm = TRUE), 5),
      hRSV_B_SD = round(sd(b, na.rm = TRUE), 6),
      F_value = round(f_val, 6),
      p_value = format(p_val, scientific = TRUE, digits = 6)
    ))
  }
}
print(results, row.names = FALSE)
# Export as CSV (saves in your working directory 
— set it if needed)
write.csv(results, "C:/Users/hp/Downloads/gene_level_comparison_ANOVA.csv", row.names = FALSE)
# ============================================================================
# GENE-LEVEL PEARSON CORRELATIONS: ENC vs GC-related metrics
# ============================================================================
# Load data (base R, no tidyverse needed)
data <- read.csv("C:/Users/hp/Downloads/COMBINED RESULTS (2).csv", stringsAsFactors = FALSE)
# Extract gene name from title column
data$Gene <- sapply(strsplit(data$title, "_"), function(x) x[1])
# Check gene names look correct
table(data$Gene)
# Metrics to correlate against ENC
gc_metrics <- c("GC", "GC1", "GC2", "GC3s", "GC12")
genes <- unique(data$Gene)
# Empty results dataframe
results <- data.frame()
for (gene in genes) {
  gene_data <- data[data$Gene == gene, ]
  for (metric in gc_metrics) {
    x <- gene_data$ENC
    y <- gene_data[[metric]]
    test <- cor.test(x, y, method = "pearson")
    results <- rbind(results, data.frame(
      Gene = gene,
      Comparison = paste0("ENC vs. ", metric),
      Pearson_r = round(test$estimate, 4),
      p_value = format(test$p.value, scientific = TRUE, digits = 4),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}
print(results, row.names = FALSE)
# Export as CSV
write.csv(results, "C:/Users/hp/Downloads/gene_level_pearson_correlations.csv", row.names = FALSE)
# ============================================================================
# GENE-LEVEL PEARSON CORRELATIONS: ENC vs GC-related metrics
# Uses existing "Gne Name" column (already extracted via Excel formula)
# ============================================================================
data <- read.csv("C:/Users/hp/Downloads/COMBINED RESULTS (2).csv", stringsAsFactors = FALSE)
# Check counts 
— should be 36 per gene (18 hRSV-A + 18 hRSV-B)
table(data$`Gne Name`)
table(data$Type)
# Metrics to correlate against ENC
gc_metrics <- c("GC", "GC1", "GC2", "GC3s", "GC12")
genes <- unique(data$`Gne Name`)
results <- data.frame()
for (gene in genes) {
  gene_data <- data[data$`Gne Name` == gene, ]
  for (metric in gc_metrics) {
    x <- gene_data$ENC
    y <- gene_data[[metric]]
    test <- cor.test(x, y, method = "pearson")
    results <- rbind(results, data.frame(
      Gene = gene,
      Comparison = paste0("ENC vs. ", metric),
      Pearson_r = round(test$estimate, 4),
      p_value = format(test$p.value, scientific = TRUE, digits = 4),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}
print(results, row.names = FALSE)
# Export as CSV
write.csv(results, "C:/Users/hp/Downloads/gene_level_pearson_correlations.csv", row.names = FALSE)
data <- read.csv("C:/Users/hp/Downloads/COMBINED RESULTS (2).csv", stringsAsFactors = FALSE)
# Check names
print(names(data))
# Check counts 
— should be 36 per gene (18 hRSV-A + 18 hRSV-B)
table(data$Gne.Name)
table(data$Type)
# Metrics to correlate against ENC
gc_metrics <- c("GC", "GC1", "GC2", "GC3s", "GC12")
genes <- unique(data$Gne.Name)
results <- data.frame()
for (gene in genes) {
  gene_data <- data[data$Gne.Name == gene, ]
  for (metric in gc_metrics) {
    x <- gene_data$ENC
    y <- gene_data[[metric]]
    test <- cor.test(x, y, method = "pearson")
    results <- rbind(results, data.frame(
      Gene = gene,
      Comparison = paste0("ENC vs. ", metric),
      Pearson_r = round(test$estimate, 4),
      p_value = format(test$p.value, scientific = TRUE, digits = 4),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}
print(results, row.names = FALSE)
# Export as CSV 
— save to Desktop instead of Downloads to avoid permission issue
write.csv(results, "C:/Users/hp/Desktop/gene_level_pearson_correlations.csv", row.names = FALSE)
# Option 2: Manually pick save location via popup
write.csv(results, file.choose(new = TRUE), row.names = FALSE)
q()
library(ggplot2)
library(ca)
# Data load
rscu_data <- read.csv("C:/CDS FASTA/RSCU_Results/RSCU_matrix_for_STATA.csv")
# RSCU matrix
rscu_matrix <- rscu_data[, 5:63]
# Correspondence Analysis
coa_result <- ca(rscu_matrix)
# Row (principal) coordinates
row_scores <- as.data.frame(cacoord(coa_result, type = "rowprincipal"))
row_scores$Sample_ID <- rscu_data$Sample_ID
row_scores$Subgroup  <- rscu_data$Subgroup
row_scores$Gene      <- rscu_data$Gene
row_scores$Country   <- rscu_data$Country
# Inertia %
inertia     <- coa_result$sv^2
inertia_pct <- round(inertia / sum(inertia) * 100, 1)
# Plot
ggplot(row_scores, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = Subgroup, shape = Gene),
             size = 2.5, alpha = 0.85) +
  stat_ellipse(aes(color = Subgroup, fill = Subgroup),
               geom = "polygon", alpha = 0.08,
               level = 0.95, type = "t") +
  scale_color_manual(values = c("hRSV-A" = "#D62728",
                                "hRSV-B" = "#1F77B4")) +
  scale_fill_manual(values  = c("hRSV-A" = "#D62728",
                                "hRSV-B" = "#1F77B4")) +
  scale_shape_manual(values = c(
    "NS1"=0, "NS2"=1, "N"=2, "P"=5, "M"=6,
    "SH"=8, "G"=15, "F"=16, "M2-1"=17, "M2-2"=18, "L"=4
  )) +
  coord_cartesian(xlim = c(-3, 4.5), ylim = c(-3, 3)) +
  labs(
    title = "Correspondence Analysis of RSCU Values in hRSV",
    x     = paste0("Dimension 1 (", inertia_pct[1], "%)"),
    y     = paste0("Dimension 2 (", inertia_pct[2], "%)"),
    color = "Subgroup", fill = "Subgroup", shape = "Gene"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
# Save
ggsave("C:/CDS FASTA/RSCU_Results/COA_Plot_final.png",
