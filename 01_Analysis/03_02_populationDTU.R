#!/usr/bin/env Rscript

# This script runs population differential transcript usage (population DTU) analysis.
# Utility functions are stored in:
# ../00_utils/00_MANTA_functions.R
#
# To identify population-associated transcript usage differences across immune lineages,
# we use a multivariate non-parametric framework to model EC ratios.
#
# Specifically, for each gene, within each immune lineage
# (CD4+ T, CD8+ T, monocytes, NK, and B cells) and each condition
# (Baseline, IAV, COV) separately, EC ratios are modeled as a function of
# technical and biological covariates (see Methods).
#
# The following model is fitted:
# 
# Ratios (EC1, ..., ECn) ~ Batch + Age + Cell Mortality + Cell type proportions + Population
#
# A significant Population term (FDR < 0.05) indicates population-associated shifts in EC relative abundances.

# Utility functions used in this script are stored in: ../00_utils/00_MANTA_functions.R


suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--rel_abundances"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'RDS object containing EC ratios (EC x Donor) ' 
  ), 
  make_option(
    c("-c", "--covariates"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'String of comma separated of covariates defining the differential model. '  # Run,Mortality,Age,POP #  Cell proportions are added from the sample metadata within this script
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file with DTU stats' 
  ), 
  make_option(
    c("-n", "--n_cores"), # 2024-09-18 Currently not useing parallel::mclapply so not much sense to include this param
    action = "store",
    default = 112,
    type = 'numeric',
    help = 'Number of cores' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## PACKAGES
suppressPackageStartupMessages(require(manta))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(reshape2))

source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_Scripts/01_DEA_common_functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")

############
## PARAMS ## 
############
rel_abundances <- opt$rel_abundances
out_file <- opt$out_file
covariates <- unlist(strsplit(opt$covariates, ","))
n_cores <- opt$n_cores - 1 # Use N -1 threads
# metadata
meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_per_Donor.ID_cell_fractions.tsv"
meta_data <- read.table(meta_data, header = TRUE, sep = "\t")

## RUN
# 1. Read input
celltype <- gsub("-.*", "", basename( rel_abundances ))
rel_abundances <- readRDS( rel_abundances )

# 2. Adjust model 
## 2.1 In case ISG activity Hallmark as covariate
hallmark.vec <- grep("Hallmark", covariates)
if(length(hallmark.vec) > 0 ){ covariates[hallmark.vec] <- paste0(covariates[hallmark.vec], "_", celltype) }
## 2.2 Add cell type fractions
cts <-  grep(paste0("^", celltype, "_"),  colnames(meta_data), value = T)
covariates <- c(covariates, cts)
covariates <- covariates[ grep("T.CD4_T.Reg", covariates, invert = T) ] # For CD4.T remove one cell type covariate because they are super correlated between them and causes MANTA to fail

# 3. MANTA
# annot_cols <- c( "gene_id", "gene_name", "transcript_name", "ec", "max.tx_name.AFB", "mean_ec_exp.AFB", "max.tx_name.EUB", "mean_ec_exp.EUB", "is.EC.switch", "MD", "tr.first", "tr.second")
manta_results_l <- manta_wrapper ( rel_abundances, meta_data, covariates )
## 3.2 Save MANTA results
saveRDS(manta_results_l, file = out_file)

# 3.3. Save empty file with model used 
parent_dir <- dirname(normalizePath(out_file))
model_file <- paste0("dtu~", paste0(gsub(" ", "", covariates), collapse = "+"), ".txt")
file.create(paste0(parent_dir,"/", model_file))