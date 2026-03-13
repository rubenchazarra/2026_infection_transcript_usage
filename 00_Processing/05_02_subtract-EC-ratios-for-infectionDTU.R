#!/usr/bin/env Rscript

# In order to conduction infection DTU analysis, subtract the relative abundances of viral-stimulated samples (e.g. IAV and COV) to baseline samples 

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--rel_abundances"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'RDS file containing relative abundance matrix for different conditions ' 
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file for difference (condition_1 - condition_2)' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

### FUNCTIONS ###
substract_matrices <- function(counts, meta_data){
  if( class(counts) == "list") { lapply( counts, function(counts) substract_matrices(counts, meta_data) ) 
    }else{
      annot_cols <- setdiff(colnames(counts), meta_data[[ "Donor.ID_Condition" ]])
      annot_df <- counts[, annot_cols]
      sample_cols <-  intersect(colnames(counts), meta_data[[ "Donor.ID_Condition" ]])
      condition.vec <-  meta_data[ match(sample_cols, meta_data[[ "Donor.ID_Condition" ]]), "Condition" ]
      condition.l <- split(sample_cols, condition.vec)
      condition.l <- lapply(condition.l, function(x) counts[, x]) # split counts matrix by Condition
      condition.l <- lapply(condition.l, function(x) setNames(x,  gsub("_.*", "", colnames(x))))
      common.donors <- Reduce(intersect, lapply(condition.l, colnames))
      condition.l <- lapply(condition.l, function(x) x[, common.donors])
      print(paste0("Substracting Condition matrices ", names(condition.l)[[1]], " - ", names(condition.l)[[2]], " ..."))
      print(paste0("N of common donors = ", length(common.donors), " ..."))
      diff_counts <- condition.l[[1]] - condition.l[[2]]
      counts <- cbind(annot_df, diff_counts)
      return(counts)
  }
}

### RUN ###

# PARAMS 
rel_abundances <- opt$rel_abundances
out_file <- opt$out_file

# 1. read
## 1.1 Metadata 
meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv"
meta_data <- read.table(meta_data, header = TRUE, sep = "\t")

## 1.2 tx2gene
# tx2gene <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv"
# tx2gene <- read.table(tx2gene, header = T, sep = "\t", na.strings = "")
# print("tx2gene file read ...")

## 1.3 Read relative abundances
celltype <- gsub("-.*", "", basename(rel_abundances))
print(paste0("Compute difference in EC abundances for ... ", celltype, " ..."))
rel_abundances <- readRDS(rel_abundances)

# 2. Substract matrices in different set ups
rel_abundances <- substract_matrices(rel_abundances, meta_data)

# 3. Save
saveRDS(rel_abundances, out_file)
print("Files saved ...")