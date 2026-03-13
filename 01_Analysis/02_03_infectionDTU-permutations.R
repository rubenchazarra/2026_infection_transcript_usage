#!/usr/bin/env Rscript

# Permutation analysis for infection-associated DTU detection.
#
# In the original (non-permuted) infection DTU analysis, for each immune lineages (CD4+ T, CD8+ T, monocytes, NK, and B cells),
# and viral response (e.g. influenza A virus (IAV) vs Baseline, SARS-CoV-2 (COV) vs Baseline),
# equivalence class (EC) differences are computed between matched infected and
# baseline samples from the same donor (e.g. Donor1_Infected − Donor1_Baseline, Donor2_Infected − Donor2_Baseline).
#
# In the permutation analysis, condition labels are randomly shuffled across
# donors, such that EC ratio differences are calculated between mismatched
# infected and baseline samples (e.g. Donor1_Infected − Donor1_Baseline, Donor2_Baseline − Donor2_Infected).
#
# This procedure is repeated 100 times, and the infection DTU analysis is run
# for each permutation using the model: 
#
# ΔEC Ratios Infected-Baseline(EC₁, ..., ECₙ) ~ Intercept + Batch + Age + Cell Mortality + Population
#
# Genes with a significant Intercept (FDR < 0.05) are considered an empirical estimate of the number of false positive infection DTU genes
#
# Utility functions used in this script are stored in: ../00_utils/00_MANTA_functions.R



suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--rel_abundances"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'RDS object containing EC ratios (not the subtraction of EC ratios between infected and baseline samples), but the actual EC ratios for infected AND baseline samples' 
  ), 
  make_option(
    c("-c", "--covariates"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'String of comma separated of covariates defining the differential model. No spaces! (e.g. Run,Mortality,Age,POP)' 
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file with DTU stats' 
  ) 
)

opt <- parse_args(OptionParser(option_list=option_list))

## FUNCTIONS 
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
## PARAMS
rel_abundances <- opt$rel_abundances
out_file <- opt$out_file
covariates <- unlist(strsplit(opt$covariates, ","))
## Metadata
meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_per_Donor.ID_cell_fractions.tsv"
meta_data <- read.table(meta_data, header = TRUE, sep = "\t")

## RUN
# 1. Read input
celltype <- gsub("-.*", "", basename(rel_abundances))
rel_abundances <- readRDS( rel_abundances )
rel_abundances <- rel_abundances[c("COV.NS", "IAV.NS" )] # 2025-08-25 Select IAV vs Baseline and COV vs Baseline, no longer interested in IAV vs COV comparison
print(paste0("Running DTU for cell type : ", celltype, " ..."))

# 2. Adjust model 
interaction <- FALSE # @RCG 2025-08-25 For infection DTU test we don't need to correct for sub-cell type fraction differences between between populations. However for the infection * population itneraction we do 
if (interaction ) {
  # 2.2 Add cell type fractions --> @
  cts <-  grep(paste0("^", celltype, "_"),  colnames(meta_data), value = T)
  covariates <- c(covariates, cts)
  ## For CD4.T remove one cell type covariate because they are super correlated between them and causes MANTA to fail. This is also done in the DGE and sQTL analysis
  covariates <- covariates[ grep("T.CD4_T.Reg", covariates, invert = T) ]
}
print(covariates)

# 3. PERMUTATIONS
iteration_l <- mantarraya_permutations_wrapper( rel_abundances , meta_data, covariates, n_permutations = 100)
## Save permutations results
saveRDS(iteration_l, out_file)

# 4. Save empty file with model used 
parent_dir <- dirname(normalizePath(out_file))
model_file <- paste0("dtu~", paste0(gsub(" ", "", covariates), collapse = "+"), ".txt")
file.create(paste0(parent_dir,"/", model_file))