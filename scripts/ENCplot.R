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
library(ggplot2)
library(dplyr)
mydata <- read.csv("C:/CDS FASTA/COMBINED RESULTS.csv", header = TRUE)
mydata$Type <- trimws(mydata$Type)
mydata$Type <- factor(mydata$Type, levels = c("hRSV-A", "hRSV-B"))
# Expected curve (Wright, 1990)
S <- seq(0.01, 0.99, by = 0.001)
curve_data <- data.frame(
  GC3s    = S,
  ENC_exp = 2 + S + 29 / (S^2 + (1 - S)^2)
)
enc_plot <- ggplot() +
  geom_line(
    data = curve_data,
    aes(x = GC3s, y = ENC_exp),
    color = "black", linewidth = 0.8
  ) +
  geom_point(
    data = mydata,
    aes(x = GC3s, y = ENC, color = Type, shape = Type),
    size = 2.5, alpha = 0.75
  ) +
  scale_color_manual(
    values = c("hRSV-A" = "#E41A1C", "hRSV-B" = "#377EB8"),
    name = "Virus Type"
  ) +
  scale_shape_manual(
    values = c("hRSV-A" = 16, "hRSV-B" = 17),
    name = "Virus Type"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(20, 65), breaks = seq(20, 65, 5)) +
  labs(
    title    = "ENC Plot Analysis of hRSV-A and hRSV-B",
    subtitle = "Points below the expected curve indicate natural selection influence",
    x = "GC3s",
    y = "ENC"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position  = "top",
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )
print(enc_plot)
ggsave("C:/CDS FASTA/ENC_Plot.png", plot = enc_plot,
