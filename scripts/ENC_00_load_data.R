# =====================================================================
# ENC_00_load_data.R
# Loads COMBINED_RESULTS.csv and assigns a gene name to every record.
# Run by the other ENC_ scripts with: source("scripts/ENC_00_load_data.R")
#
# Input : data/COMBINED_RESULTS.csv   (396 records: 18 strains x 11 genes x 2 subgroups)
# =====================================================================

dir.create("results", showWarnings = FALSE)

dat <- read.csv("data/COMBINED_RESULTS.csv", stringsAsFactors = FALSE)

# Gene name from the sequence title (checked in this order)
x <- dat$title
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

stopifnot(!any(dat$Gene == "UNKNOWN"))   # every record must be assigned
print(table(dat$Gene, dat$Type))         # expect 18 per gene per subgroup

# Gene order used in the tables
genes <- c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2-1", "M2-2", "L")
