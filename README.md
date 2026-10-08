# Codon-Usage-Bias
# Codon-Usage-Bias

R scripts and data for the codon usage bias analysis of hRSV-A and hRSV-B
(thesis, FJWU Rawalpindi).

## Folders
- DATA/    : input files (COMBINED_RESULTS.csv, RSCU_matrix_for_STATA.csv,
             CAI_HRSV_A_and_B_Combined.xlsx, accession numbers)
- scripts/ : R scripts
- results/ : created automatically when scripts are run

## How to run
Set the R working directory to the project folder (the one containing
DATA/ and scripts/), then run the scripts in this order:
1. ENC_01_descriptive.R  (Tables 3.1, 3.5)
2. ENC_02_anova.R        (Table 3.2)
3. ENC_03_pearson.R      (Tables 3.3, 3.4)
4. ENCplot.R             (ENC plot)
5. COA_plot.R            (correspondence analysis)
6. RSCU_table.R          (see note below)

## Note
RSCU_table.R needs CDS FASTA files, which are not included. Download the
sequences from NCBI using the accession numbers in DATA/, extract the CDS
of each gene, and place them in DATA/CDS_FASTA/hRSV-A/ and hRSV-B/.

## Requirements
R, with packages: tidyverse, ggplot2, ca, readxl.
