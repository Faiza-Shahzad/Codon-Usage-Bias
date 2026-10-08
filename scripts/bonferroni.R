# =====================================================================
# 03_bonferroni.R
# Bonferroni correction for the 11 gene-wise t-tests.
# Corrected p = raw p x 11 (capped at 1). Significant if corrected p < 0.05.
# Output: results/CAI_bonferroni.csv
# =====================================================================

source("scripts/00_load_data.R")

pv <- do.call(rbind, lapply(genes, function(g) {
  a <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-A"]
  b <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-B"]
  data.frame(Gene = g,
             p_student = t.test(a, b, var.equal = TRUE)$p.value,
             p_welch   = t.test(a, b)$p.value)
}))

pv$bonf_student <- p.adjust(pv$p_student, method = "bonferroni")
pv$bonf_welch   <- p.adjust(pv$p_welch,   method = "bonferroni")
pv$significant_student <- ifelse(pv$bonf_student < 0.05, "Yes", "No")
pv$significant_welch   <- ifelse(pv$bonf_welch   < 0.05, "Yes", "No")

print(format(pv, digits = 4), row.names = FALSE)
write.csv(pv, "results/CAI_bonferroni.csv", row.names = FALSE)
