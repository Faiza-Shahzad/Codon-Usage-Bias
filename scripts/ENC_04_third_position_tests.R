# =====================================================================
# ENC_04_third_position_tests.R
# Table 3.6: comparison of third-position nucleotide frequencies
# (A3s, T3s, C3s, G3s) between hRSV-A and hRSV-B across all 396 records
# (198 per subgroup): t-test and Mann-Whitney U test.
#
# t.test() in R uses Welch's test by default; this reproduces the
# reported t and p values. t is hRSV-A minus hRSV-B.
# Output: results/Table3_6_third_position_tests.csv
# =====================================================================

source("scripts/ENC_00_load_data.R")

vars <- c("A3s", "T3s", "C3s", "G3s")

t3_6 <- do.call(rbind, lapply(vars, function(v) {
  a  <- dat[[v]][dat$Type == "hRSV-A"]
  b  <- dat[[v]][dat$Type == "hRSV-B"]
  tt <- t.test(a, b)                       # Welch's t-test
  mw <- wilcox.test(a, b, exact = FALSE)   # Mann-Whitney U (ties present)
  data.frame(Parameter = v,
             t_value = unname(tt$statistic),
             p_t_test = tt$p.value,
             significant_t_test = ifelse(tt$p.value < 0.05, "Yes", "No"),
             p_MannWhitney = mw$p.value,
             significant_MannWhitney = ifelse(mw$p.value < 0.05, "Yes", "No"))
}))

print(format(t3_6, digits = 4), row.names = FALSE)
write.csv(t3_6, "results/Table3_6_third_position_tests.csv", row.names = FALSE)
