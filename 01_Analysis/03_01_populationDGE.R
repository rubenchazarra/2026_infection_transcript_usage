#!/usr/bin/env Rscript 

# This script performs population differential gene expression (population DEG) analysis
# between individuals from different populations (e.g. European vs African).
#
# For each immune lineage (CD4+ T, CD8+ T, monocytes, NK, and B cells), in each separate condition (Baseline, IAV, COV)
# gene expression is modeled as a function of technical and biological covariates
# using a linear mixed-effects model.
#
# The following model is fitted for each gene:
#   Expression ~ (1 | Batch) + (1 | Library) + Age + Cell Mortality + Cell type proportions + Population
#
# The Population term captures gene expression differences between populations.
# Genes with a significant Population coefficient (FDR < 0.01) and | log2FC | > 0.2 are considered population-associated DEGs.

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'CPM Normalized Expression Counts' 
  ), 
  make_option(
    c("-c", "--covariates"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'String of comma separated of covariates defining the differential model. No spaces!' 
  ), 
  make_option(
    c("-n", "--n_cores"),
    action = "store",
    default = 112,
    type = 'numeric',
    help = 'Number of cores' 
  ),
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## FUNCTIONS ##
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_Scripts/01_DEA_common_functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/01_DGE/00_Scripts/00_DGE-functions.R")

###############
##### RUN ##### 
###############

## PARAMS
counts <- opt$counts
covariates <- unlist(strsplit(opt$covariates, ","))
n_cores <- opt$n_cores
out_file <- opt$out_file

# PARALLELIZATION
param = SnowParam(n_cores, "SOCK", progressbar=TRUE)
register(param)

# 1. Read counts
celltype <- gsub("-.*", "", basename(counts))
counts <- readRDS(counts)

# 2. METADATA
sample_meta_data <- load.Aquino.metadata()

# 3. MODEL
## 3.1 Add cell type fractions as covaraites
cts <-  grep(paste0("^", celltype, "_"),  colnames(sample_meta_data), value = T)
covariates <- c(covariates, cts)
#  # If we don't do this, we lead to the following error logs/24110056_4-popDGE.err. WE did the same for sQTLsseekeR2
if(celltype %in% "T.CD4"){  covariates <- grep("T.CD4_T.Reg", covariates, value = T, invert = T) } 

# 4. DGE
dge.l <- dge.limma.dream.wrapper( counts, sample_col = "Donor.ID_Condition", sample_meta_data, covariates, param )

# 5. Save output
saveRDS(dge.l, out_file)
print("Output file saved ...")

# 5.2 Save empty file with model
parent_dir <- dirname(out_file)
model_file <- paste0("dge~", paste0(gsub(" ", "", covariates), collapse = "+"), ".txt")
file.create(paste0(parent_dir,"/", model_file))
