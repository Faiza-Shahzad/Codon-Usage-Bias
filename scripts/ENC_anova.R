# =====================================================================
# ENC_02_anova.R
# Table 3.2: one-way ANOVA of ENC between hRSV-A and hRSV-B, per gene
# (18 strains per subgroup, so n = 36 per gene).
# Output: results/Table3_2_anova_ENC.csv
# =====================================================================

source("scripts/ENC_00_load_data.R")

anova_tab <- do.call(rbind, lapply(genes, function(g) {
  d   <- dat[dat$Gene == g, ]
  fit <- summary(aov(ENC ~ Type, data = d))[[1]]
  a   <- d$ENC[d$Type == "hRSV-A"]
  b   <- d$ENC[d$Type == "hRSV-B"]
  data.frame(Gene = g,
             A_mean = mean(a), A_SD = sd(a),
             B_mean = mean(b), B_SD = sd(b),
             F_value = fit[["F value"]][1],
             p_value = fit[["Pr(>F)"]][1])
}))
anova_tab$significant <- ifelse(anova_tab$p_value < 0.05, "Yes", "No")

print(format(anova_tab, digits = 4), row.names = FALSE)
write.csv(anova_tab, "results/Table3_2_anova_ENC.csv", row.names = FALSE)
