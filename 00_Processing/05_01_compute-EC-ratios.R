#!/usr/bin/env Rscript

# Compute equivalence class (EC) relative abundances for later modelling of differential transcript usage (DTU) differences upon infection, between populations, or infection*population interactions

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-c", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'RDS or TXT object with EC x cell matrix' 
  ), 
  make_option(
    c("-t", "--min.transcript.exp"),
    action = "store",
    default = 1,
    type = 'numeric',
    help = 'sQTLseekeR2::prepare.trans.exp( min.transcript.exp = min.transcript.exp)' 
  ), 
  make_option(
    c("-g", "--min.gene.exp"),
    action = "store",
    default = 5,
    type = 'numeric',
    help = 'sQTLseekeR2::prepare.trans.exp( min.gene.exp = min.gene.exp)' 
  ), 
  make_option(
    c("-p", "--min.prop"),
    action = "store",
    default = 0.4,
    type = 'numeric',
    help = 'sQTLseekeR2::prepare.trans.exp( min.prop = min.prop)' 
  ), 
  make_option(
    c("-d", "--min.dispersion"),
    action = "store",
    default = 0.01,
    type = 'numeric',
    help = 'sQTLseekeR2::prepare.trans.exp( min.dispersion = min.dispersion)' 
  ), 
  make_option(
    c("-m", "--meta_col"),
    action = "store",
    default = "POP",
    type = 'character',
    help = 'Metadata column to compute MD (DTU effect size) by: one of "Condition" (for infection DTU) or "POP" (for population DTU)' 
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file with EC Relative Abundances'
  ) 
)

opt <- parse_args(OptionParser(option_list=option_list))

### FUNCTIONS ###
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_Scripts/01_DEA_common_functions.R")

# PARAMS 
counts <- opt$counts
min.transcript.exp <- opt$min.transcript.exp
min.gene.exp <- opt$min.gene.exp
min.prop <- opt$min.prop
min.dispersion <- opt$min.dispersion
out_file <- opt$out_file
n_cores <- 56
meta_col <- opt$meta_col

# 0. Print param combination 
print(paste0("Parameter combination: --min.transcript.exp=",min.transcript.exp, "  --min.gene.exp=", min.gene.exp, "  --min.prop=", min.prop, "  --min.dispersion=", min.dispersion))
# 1. Data 
## 1.1 Metadata
meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv"
meta_data <- read.table(meta_data, header = TRUE, sep = "\t")

## 1.2 tx2gene
tx2gene <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv"
tx2gene <- read.table(tx2gene, header = T, sep = "\t", na.strings = "")
print("tx2gene file read ...")

# 1.3 Read TCC counts
celltype <- gsub("-.*", "", basename(counts))
counts <- readRDS(counts)
print(paste0("Compute EC rel abundances for: ", celltype, " ..."))

# 2. Relative abundances
rel.counts <- sQTLseekeR2.rel_abundances.wrapper(counts,
                                                 tx2gene, 
                                                 meta_data, 
                                                 n_cores, 
                                                 min.transcript.exp=min.transcript.exp, 
                                                 min.gene.exp=min.gene.exp, 
                                                 min.prop=min.prop, 
                                                 min.dispersion= min.dispersion, 
                                                 meta_col = meta_col
                                                 )
print("Relative abundances computed ...")

# 3. Save list output
saveRDS(rel.counts, out_file)
print("File saved ...")
