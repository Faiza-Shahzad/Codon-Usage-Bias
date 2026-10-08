# =====================================================================
# 00_load_data.R
# Loads the CAI data and assigns a gene name to every sequence.
# Run by the other scripts with: source("scripts/00_load_data.R")
#
# Input : data/CAI_HRSV_A_and_B_Combined.xlsx  (sheet "CAI_Results")
# Requires: readxl  ->  install.packages("readxl")
# =====================================================================

library(readxl)
dir.create("results", showWarnings = FALSE)

dat <- as.data.frame(read_excel("DATA/CAI_HRSV_A_and_B_Combined.xlsx",
                                sheet = "CAI_Results"))
x <- dat$Genes
dat$Gene <- "UNKNOWN"
dat$Gene[startsWith(x, "L_")]  <- "L"
dat$Gene[startsWith(x, "F_")]  <- "F"
dat$Gene[startsWith(x, "G_")]  <- "G"
dat$Gene[startsWith(x, "SH")]  <- "SH"
dat$Gene[startsWith(x, "M_")]  <- "M"
dat$Gene[startsWith(x, "P_")]  <- "P"
dat$Gene[startsWith(x, "N_")]  <- "N"
dat$Gene[startsWith(x, "NS2")] <- "NS2"
dat$Gene[startsWith(x, "NS1")] <- "NS1"
dat$Gene[grepl("M2-2", x)]     <- "M2-2"
dat$Gene[grepl("M2-1", x)]     <- "M2-1"

stopifnot(!any(dat$Gene == "UNKNOWN"))   
print(table(dat$Gene, dat$Subgroup))     

genes <- c("F", "G", "L", "M", "M2-1", "M2-2", "N", "NS1", "NS2", "P", "SH")
