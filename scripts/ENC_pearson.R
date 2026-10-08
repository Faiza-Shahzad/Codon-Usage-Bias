# =====================================================================
# ENC_03_pearson.R
# Tables 3.3 (hRSV-A) and 3.4 (hRSV-B): Pearson correlation between ENC
# and each GC parameter, per gene, within each subgroup (n = 18 strains).
#
# Bonferroni: p-values are corrected within each gene and subgroup, over
# the comparisons that could be computed (5, or 4 when one GC parameter
# is constant across all strains, which is reported as "Invariant").
# Output: results/Table3_3_3_4_pearson.csv
# =====================================================================

source("scripts/ENC_00_load_data.R")

gc_vars <- c("GC", "GC1", "GC2", "GC3s", "GC12")

pearson_tab <- do.call(rbind, lapply(c("hRSV-A", "hRSV-B"), function(tp) {
  do.call(rbind, lapply(genes, function(g) {
    d <- dat[dat$Gene == g & dat$Type == tp, ]

    out <- do.call(rbind, lapply(gc_vars, function(v) {
      if (sd(d[[v]]) == 0 || sd(d$ENC) == 0) {       # constant variable: r undefined
        return(data.frame(Subgroup = tp, Gene = g,
                          Comparison = paste("ENC vs.", v),
                          Pearson_r = NA, p_value = NA))
      }
      ct <- cor.test(d$ENC, d[[v]], method = "pearson")
      data.frame(Subgroup = tp, Gene = g,
                 Comparison = paste("ENC vs.", v),
                 Pearson_r = unname(ct$estimate), p_value = ct$p.value)
    }))

    n_valid <- sum(!is.na(out$p_value))
    out$p_bonferroni <- p.adjust(out$p_value, method = "bonferroni", n = n_valid)
    out$significant  <- ifelse(is.na(out$p_value), "Invariant",
                        ifelse(out$p_bonferroni < 0.05, "Yes", "No"))
    out
  }))
}))

print(format(pearson_tab, digits = 4), row.names = FALSE)
write.csv(pearson_tab, "results/Table3_3_3_4_pearson.csv", row.names = FALSE)
