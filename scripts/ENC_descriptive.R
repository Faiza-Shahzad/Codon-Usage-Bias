# =====================================================================
# ENC_01_descriptive.R
# Mean and SD by subgroup.
#   Table 3.1: ENC and GC-content parameters (GC values as percentages)
#   Table 3.5: third-position nucleotides and GC parameters (as fractions)
# Outputs: results/Table3_1_descriptive.csv, results/Table3_5_third_position.csv
# =====================================================================

source("scripts/ENC_00_load_data.R")

describe <- function(vars, scale = 1) {
  do.call(rbind, lapply(vars, function(v) {
    a <- dat[[v]][dat$Type == "hRSV-A"] * scale
    b <- dat[[v]][dat$Type == "hRSV-B"] * scale
    data.frame(Parameter = v,
               hRSV_A_mean = mean(a), hRSV_A_SD = sd(a),
               hRSV_B_mean = mean(b), hRSV_B_SD = sd(b))
  }))
}

# Table 3.1: ENC as is, GC parameters x 100 (percent)
t3_1 <- rbind(describe("ENC"),
              describe(c("GC", "GC1", "GC2", "GC3s", "GC12"), scale = 100))

# Table 3.5: fractions
t3_5 <- describe(c("A3s", "T3s", "C3s", "G3s", "GC", "GC1", "GC2", "GC3s", "GC12"))

cat("\nTable 3.1\n"); print(format(t3_1, digits = 4), row.names = FALSE)
cat("\nTable 3.5\n"); print(format(t3_5, digits = 4), row.names = FALSE)

write.csv(t3_1, "results/Table3_1_descriptive.csv",   row.names = FALSE)
write.csv(t3_5, "results/Table3_5_third_position.csv", row.names = FALSE)
