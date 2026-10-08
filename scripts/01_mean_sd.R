# =====================================================================
# 01_mean_sd.R
# Mean and standard deviation of CAI per gene and subgroup,
# plus overall means across the 11 genes.
# Output: results/CAI_mean_sd.csv
# =====================================================================

source("scripts/00_load_data.R")

desc <- do.call(rbind, lapply(genes, function(g) {
  a <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-A"]
  b <- dat$CAI[dat$Gene == g & dat$Subgroup == "hRSV-B"]
  data.frame(Gene = g,
             n_A = length(a), n_B = length(b),
             mean_A = mean(a), mean_B = mean(b),
             sd_A = sd(a), sd_B = sd(b),
             diff_A_B = mean(a) - mean(b))
}))

print(format(desc, digits = 4), row.names = FALSE)

cat("\nMean CAI across all genes\n")
cat("hRSV-A:", round(mean(desc$mean_A), 4), "\n")
cat("hRSV-B:", round(mean(desc$mean_B), 4), "\n")
cat("Difference:", round(mean(desc$mean_A) - mean(desc$mean_B), 4), "\n")

write.csv(desc, "results/CAI_mean_sd.csv", row.names = FALSE)
