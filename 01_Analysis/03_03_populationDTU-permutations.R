#!/usr/bin/env Rscript


# Permutation analysis for population-associated DTU detection.
#
# This script evaluates the false positive rate of the population DTU framework
# by randomly shuffling donor labels within each immune lineage and condition
# and repeating the DTU analysis 100 times.
#
# For each immune lineage (e.g. CD4+ T cells) and each condition
# (Baseline, IAV, COV):
# 1) randomly permute donor labels across samples;
# 2) re-run the population DTU analysis using the permuted labels;
# for each gene we run the following model Ratios (EC1, ..., ECn) ~ Batch + Age + Cell Mortality + Cell type proportions + Population
# 3) record the number of genes with a significant Population term
#    (FDR < 0.05, Benjamini–Hochberg correction).
#
# This analysis provides an empirical estimate of the number of false positive
# population DTU genes expected under the null hypothesis of no association
# between population labels and EC ratio variation.
#
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
    help ='String of comma separated of covariates defining the differential model. '  # Run,Mortality,Age,POP #  Cell proportions are added from the sample metadata within this script
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file with DTU results' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## PACKAGES
suppressPackageStartupMessages(require(manta))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(reshape2))

## FUNCTIONS
# source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_Scripts/01_DEA_common_functions.R")
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
celltype <- gsub("-.*", "", basename( rel_abundances ))
rel_abundances <- readRDS( rel_abundances )
print(paste0("Running DTU permutations for cell type : ", celltype, " ..."))

# 2. Adjust model 
## 2.1 For each celltype, add sub-cell type abundances as covariates in the model. This is done to correct for population cell composition differences, which are reported in Aquino et al., 2023
cts <-  grep(paste0("^", celltype, "_"),  colnames(meta_data), value = T)
covariates <- c(covariates, cts)
covariates <- covariates[ grep("T.CD4_T.Reg", covariates, invert = T) ] # For CD4.T remove one sub-cell type covariate because they are super correlated between them and causes MANTA to fail
model <- as.formula(paste("dtu ~ ", paste(c(covariates), collapse = " + ")))
print( model )

# 3. Run permutations
iteration_l <- manta_permutations_wrapper(rel_abundances, meta_data, covariates, n_permutations = 100)
  
## 3.2  Save permutations results
saveRDS(iteration_l, out_file)
