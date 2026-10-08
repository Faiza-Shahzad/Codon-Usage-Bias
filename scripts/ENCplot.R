library(tidyverse)
library(ggplot2)
library(dplyr)

dir.create("results", showWarnings = FALSE)

# Load data
mydata <- read.csv("data/COMBINED_RESULTS.csv", header = TRUE)
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
ggsave("results/ENC_Plot.png", plot = enc_plot, width = 8, height = 6, dpi = 300)
