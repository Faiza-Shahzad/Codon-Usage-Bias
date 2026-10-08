# =====================================================================
# RSCU_table.R
# Relative Synonymous Codon Usage (RSCU) for hRSV-A and hRSV-B,
# pooled over all coding sequences (CDS) of each subgroup.
#
# RSCU = observed count of a codon /
#        (total count of its amino acid / number of synonymous codons)
# RSCU = 1 means no preference; > 1 preferred; < 1 avoided.
#
# Input : FASTA files with the CDS of each gene, one folder per subgroup:
#           DATA/CDS_FASTA/hRSV-A/*.fasta
#           DATA/CDS_FASTA/hRSV-B/*.fasta
#         (FASTA files are not included; download the sequences from NCBI
#          using the accession numbers provided in DATA/)
# Output: results/RSCU_final_table.csv
#
# Base R only, no extra packages needed.
# AUG (Met), UGG (Trp) and stop codons are excluded (59 codons remain).
# =====================================================================

fasta_dir  <- "DATA/CDS_FASTA"
subgroups  <- c("hRSV-A", "hRSV-B")
dir.create("results", showWarnings = FALSE)

# ---- 1. Standard genetic code ---------------------------------------
bases <- c("T", "C", "A", "G")
g <- expand.grid(third = bases, second = bases, first = bases,
                 stringsAsFactors = FALSE)
all_codons <- paste0(g$first, g$second, g$third)
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]

aa3 <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", Q = "Gln",
         E = "Glu", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
         F = "Phe", P = "Pro", S = "Ser", T = "Thr", Y = "Tyr", V = "Val")

# ---- 2. Helper functions --------------------------------------------
read_fasta <- function(file) {
  lines <- trimws(readLines(file, warn = FALSE))
  lines <- lines[lines != ""]
  is_header <- startsWith(lines, ">")
  id <- cumsum(is_header)
  seqs <- tapply(lines[!is_header], id[!is_header], paste, collapse = "")
  toupper(gsub("U", "T", unname(seqs), ignore.case = TRUE))
}

count_codons <- function(seqs) {
  counts <- setNames(numeric(64), all_codons)
  for (s in seqs) {
    n <- nchar(s) %/% 3
    if (n == 0) next
    codons <- substring(s, seq(1, by = 3, length.out = n),
                           seq(3, by = 3, length.out = n))
    tab <- table(codons[codons %in% all_codons])
    counts[names(tab)] <- counts[names(tab)] + as.numeric(tab)
  }
  counts
}

calc_rscu <- function(counts) {
  tab <- data.frame(Codon = all_codons, AA = aa,
                    Count = as.numeric(counts[all_codons]),
                    stringsAsFactors = FALSE)
  tab <- tab[!(tab$AA %in% c("M", "W", "*")), ]          # keep 59 codons
  tab <- tab[order(match(tab$AA, unique(tab$AA))), ]       # group by amino acid
  n_syn <- ave(tab$Count, tab$AA, FUN = length)
  total <- ave(tab$Count, tab$AA, FUN = sum)
  tab$RSCU <- tab$Count / (total / n_syn)
  tab
}

# ---- 3. Compute RSCU for each subgroup ------------------------------
res <- list()
for (sg in subgroups) {
  files <- list.files(file.path(fasta_dir, sg),
                      pattern = "\\.(fa|fasta|fna|txt)$",
                      full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(files) > 0)
  seqs <- unlist(lapply(files, read_fasta))
  cat(sg, ":", length(seqs), "sequences read\n")      # expect 198 per subgroup
  res[[sg]] <- calc_rscu(count_codons(seqs))
}

final <- data.frame(
  `Amino Acid` = aa3[res[[1]]$AA],
  Codon        = gsub("T", "U", res[[1]]$Codon),
  `hRSV-A RSCU` = round(res[["hRSV-A"]]$RSCU, 3),
  `hRSV-B RSCU` = round(res[["hRSV-B"]]$RSCU, 3),
  check.names = FALSE, row.names = NULL
)

# ---- 4. Check: RSCU values of an amino acid must sum to its codon number
for (sg in subgroups) {
  s <- tapply(res[[sg]]$RSCU, res[[sg]]$AA, sum)
  n <- as.numeric(table(res[[sg]]$AA)[names(s)])
  stopifnot(all(abs(s - n) < 1e-6))
}

# ---- 5. Optional: compare with your existing table ------------------
if (file.exists("DATA/RSCU_final_table.csv")) {
  old <- read.csv("DATA/RSCU_final_table.csv", check.names = FALSE)
  old <- old[match(final$Codon, old$Codon), ]
  cat("Max difference vs DATA/RSCU_final_table.csv\n")
  cat(" hRSV-A:", max(abs(old[["hRSV-A RSCU"]] - final[["hRSV-A RSCU"]])), "\n")
  cat(" hRSV-B:", max(abs(old[["hRSV-B RSCU"]] - final[["hRSV-B RSCU"]])), "\n")
}

print(final, row.names = FALSE)
write.csv(final, "results/RSCU_final_table.csv", row.names = FALSE)
