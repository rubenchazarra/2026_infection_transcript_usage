#!/usr/bin/env Rscript 

# This script splits each pseudo-bulked and filtered immune lineage object
# (e.g. CD4+ T, CD8+ T, monocytes, NK, and B cells) by Condition
# (Baseline, IAV, COV) and generates two sets of objects:
#
# 1) split_Condition suffix:
#    Each immune lineage is saved as a list of individual conditions
#    (Baseline, IAV, COV). These objects are used to run population-level
#    differential transcript usage (DTU), differential gene expression (DGE),
#    and transcript usage QTL (TU-QTL) analyses separately for each
#    lineage–condition combination (e.g. CD4+T-Baseline).
#
# 2) group_Condition suffix:
#    For each immune lineage, viral responses are grouped into contrasts
#    (IAV vs Baseline, COV vs Baseline, and IAV vs COV). These objects are used
#    to perform infection-associated DTU, DGE, and infection × population
#    interaction analyses.

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Raw Expression Counts' 
  ), 
  make_option(
    c("-o", "--out_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to save output file' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

##### FUNCTIONS ##### 

split_by <- function(counts, meta_data, split_column = "Condition"){
  # Split counts dataframe by split_column in meta_data
  ## checks
  if(!split_column %in% names(meta_data)) stop (paste0(split_col, " ---- Not in meta-data file"))
  split_levels <- levels(factor(meta_data [[ split_column ]] ))
  ## discriminate column identifying sample: Sample.ID or Donor.ID_Condition ## 2024-08-21 This differentiate if the input files have been averaged across replicates ( Donor.ID_Condition ) or not ( Sample.ID )
  sample_col <- "Donor.ID_Condition" # Data comes from averaged replicates
  
  # Iterate
  if (class(counts) == "list"){ 
    lapply(counts, function(df) split_by(df, meta_data, split_column) )
    
  }else{
    ## split
    counts.l <- lapply(split_levels, function(x) { 
      counts[ , colnames(counts) %in% meta_data[meta_data[[split_column]] == x, sample_col ] ] 
    } )
    names(counts.l) <- split_levels
    ## remove empty levels
    counts.l <- counts.l[unlist(lapply(counts.l, function(x) ncol(x) > 0))]
    return(counts.l)
  }
}

group_by_Condition <- function(counts.l){
  
  # Group list of Condition count dataframes ready for pairwise comparison: "COV-NS", "IAV-NS", "IAV-COV"
  groups <- list("COV.NS" = c("COV", "NS"), "IAV.NS" = c("IAV", "NS"), "IAV.COV" = c("IAV", "COV"))
  
  if(class(counts.l[[1]]) == "list"){ # if object inside list is list, iterate. This is for when we have data split by POP
    lapply(counts.l, group_by_Condition)
  }else{
    counts.l <- lapply(groups, function(x) {
      l <- list(counts.l[[ x[1] ]], counts.l[[ x[2] ]]) 
      setNames(l, nm = x)
    })
    names(counts.l) <- names(groups)
    # cbind to one dataframe per pairwise comparison --> IMPORTANT TO have the data merged for CPM feature filtering 
    counts.l <- lapply(counts.l, function(l) Reduce(cbind, l ) )
    return(counts.l)
  }
}

##### RUN ##### 

## PARAMS
counts <- opt$counts
out_path <- opt$out_path
out_file <- paste0(out_path, "/", basename(counts))
## meta data
meta_data <- "../../00_MetaData/sample_metadata_expanded.tsv"
meta_data <- read.table(meta_data, header = TRUE, sep = "\t")

# 0. Read input file 
counts <- readRDS(counts)
print("File read ...")

# 1. Split by Condition 
counts.l <- split_by( counts, meta_data = meta_data, split_column = "Condition" )
## 1.2 Save
out_file <- gsub( ".rds", "-split-Condition.rds", out_file )
saveRDS( counts.l, out_file )

# 2. Group by Condition
counts.l <- group_by_Condition( counts.l )
## 2.2 Save
out_file <- gsub( "-split-Condition.rds", "-group-Condition.rds", out_file )
saveRDS( counts.l, out_file )
print("File saved ...")
