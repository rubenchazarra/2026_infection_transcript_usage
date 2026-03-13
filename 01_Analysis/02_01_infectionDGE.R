#!/usr/bin/env Rscript 

# This script performs population differential gene expression (population DEG) analysis
# between samples from different viral stimullations (e.g. IAV vs Baseline, COV vs Baseline ).
#
# For each immune lineage (CD4+ T, CD8+ T, monocytes, NK, and B cells), in each viral response (IAV vs Baseline, COV vs Baseline)
# gene expression is modeled as a function of technical and biological covariates
# using a linear mixed-effects model.
#
# The following model is fitted for each gene:
#   Expression ~ (1| Batch) + (1|Library) +  Age +  Cell Mortality + Population + Condition +  (1|Donor). 
#
# The Population term captures gene expression differences between populations.
# Genes with a significant Condition coefficient (FDR < 0.01) and | log2FC | > 0.5 are considered
# infection-associated DEGs.

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'CPM Normalized Pseudo-bulked Gene Expression Counts' 
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

# ADD individual_ID as random effect
if("Run" %in% covariates) { print("Adding 'Run' as random effect ..."); covariates[ grep("Run", covariates) ] <- "(1|Run)" }
if("Library" %in% covariates) { print("Adding 'Library' as random effect ..."); covariates[ grep("Library", covariates) ] <- "(1|Library )" }
if("Donor.ID" %in% covariates) { print("Adding 'Donor.ID' as random effect ..."); covariates[ grep("Donor.ID", covariates) ] <- "(1|Donor.ID)" }

# 1. Read counts
celltype <- gsub("-.*", "", basename(counts))
counts <- readRDS(counts)

# 2. DGE MODEL
model <- as.formula(paste(" ~  ", paste(c(covariates), collapse = " + ")))
print(paste0("Running model: ", paste0(model[1], model[2])))

# 3. METADATA
sample_meta_data <- load.Aquino.metadata()

# 4. DGE
dge.l <- dge.limma.dream.wrapper(counts, sample_meta_data, model, param)

# 5. Save output
saveRDS(dge.l, out_file)
print("Output file saved ...")

# 5.2 Save empty file with model
parent_dir <- dirname(out_file)
model_file <- paste0("dge~", paste0(gsub(" ", "", covariates), collapse = "+"), ".txt")
file.create(paste0(parent_dir,"/", model_file))
