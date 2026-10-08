library(ggplot2)
library(ca)          # install.packages("ca") if needed

dir.create("results", showWarnings = FALSE)

# Data load
rscu_data <- read.csv("DATA/RSCU_matrix_for_STATA.csv")

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
coa_plot <- ggplot(row_scores, aes(x = Dim1, y = Dim2)) +
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

print(coa_plot)

# Save
ggsave("results/COA_Plot_final.png", plot = coa_plot,
       width = 9, height = 6, dpi = 300)
