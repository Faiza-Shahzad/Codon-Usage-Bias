# =====================================================================
# 02_t_tests.R
# t-test of CAI between hRSV-A and hRSV-B for each gene.
# Student's t-test (equal variances) and Welch's t-test are both given.
# Use ONE of them throughout and state which in your methods.
# Output: results/CAI_t_tests.csv
# =====================================================================

source("scripts/00_load_data.R")

tt <- do.call(rbind, lapply(genes, function(g) {
  a <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-A"]
  b <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-B"]
  student <- t.test(a, b, var.equal = TRUE)   # Student's t-test
  welch   <- t.test(a, b)                     # Welch's t-test (R default)
  data.frame(Gene = g,
             t_student = unname(student$statistic), p_student = student$p.value,
             t_welch   = unname(welch$statistic),   p_welch   = welch$p.value)
}))

print(format(tt, digits = 4), row.names = FALSE)
write.csv(tt, "results/CAI_t_tests.csv", row.names = FALSE)
