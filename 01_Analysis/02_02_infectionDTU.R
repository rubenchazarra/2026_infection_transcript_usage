#!/usr/bin/env Rscript

# This script runs infection differential transcript usage (infection DTU) analysis. 
#
# To identify infection-associated transcript usage differences across immune lineages (CD4+ T, CD8+ T, monocytes, NK, and B cells),
# and viral responses (IAV vs Baseline, COV vs Baseline)
# we developed a strategy based on multivariate non-parametric modelling of EC ratios.
# Specifically, for each donor we subtracted EC relative abundances between infected
# and baseline samples (e.g. Donor1_IAV − Donor1_Baseline) and jointly modelled the
# resulting EC difference vectors for each gene (see Methods).
#
# This approach captures coordinated shifts in transcript usage while controlling
# for technical and biological covariates (Fig. 2a).
#
# The following model is applied for each gene:
#   ΔEC Ratios Infected-Baseline(EC₁, ..., ECₙ) ~ Intercept + Batch + Age + Cell Mortality + Population
#
# The Intercept term represents the infection effect (infected − baseline).
# The remaining covariates capture interactions between infection and each
# technical or biological factor included in the model.
#
# Genes with a significant Intercept (FDR < 0.05) are considered infection-DTU
#
# Functions are stored in the ../00_utils/00_MANTA_functions.R script.

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--rel_abundances"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'RDS object containing substraction of EC rations between infected vs baseline samples' 
  ), 
  make_option(
    c("-c", "--covariates"),
    action = "store",
    default = "",
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

## PACKAGES
suppressPackageStartupMessages(require(manta))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(reshape2))

## FUNCTIONS
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")

## PARAMS
rel_abundances <- opt$rel_abundances
out_file <- opt$out_file
covariates <- unlist(strsplit(opt$covariates, ","))
## metadata
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

# 2. Run MANTARRAYA wrapper
manta_results_l <- mantarraya_meta_wrapper ( rel_abundances, meta_data, covariates )
## 2.2 Save MANTA results
saveRDS(manta_results_l, file = out_file)
# 2.3 Save empty file with model used 
parent_dir <- dirname(normalizePath(out_file))
model_file <- paste0("dtu~", paste0(gsub(" ", "", covariates), collapse = "+"), ".txt")
file.create(paste0(parent_dir,"/", model_file))
