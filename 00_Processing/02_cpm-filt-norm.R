#!/usr/bin/env Rscript 

# This script filters equivalence classes (ECs) or genes from the pseudo-bulk object based on CPM normalization 
# Specifically we retain all ECs or genes expressed at >= 1 mean CPM across all samples

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
    c("-m", "--cpm.filt"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = 'CPM filtering' 
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

##### FUNCTIONS ##### 

## Filter by CPM 
# filter_by_cpm_counts <-  function(counts, cpm.filt = 1, frac.sample.filt  = 0.4){
#   suppressPackageStartupMessages(require(edgeR))
#   if(class(counts) == "list"){lapply(counts, function(counts) filter_by_cpm_counts( counts, cpm.filt , frac.sample.filt ) )
#     }else{
#   # Filter by CPM expression
#   n.sample.filt <- round(frac.sample.filt * ncol(counts))
#   print(paste0("Filtering features with < ", cpm.filt, " CPMs in < ", n.sample.filt, " samples ..."))
#   ## features to retain
#   isexpr = rowSums( cpm(counts, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0) > cpm.filt) >= n.sample.filt
#   print(paste0("N features expressed: ", table(isexpr)["TRUE"], " / ", length(isexpr), " ..."))
#   counts[isexpr,]
#     }
# }

filter_by_mean_cpm_counts <-  function(counts, cpm.filt = 1){
  suppressPackageStartupMessages(require(edgeR))
  if(class(counts) == "list"){lapply(counts, function(counts) filter_by_mean_cpm_counts( counts, cpm.filt ) )
  }else{
    # Filter by Mean CPM expression
    print(paste0("Filtering features with mean CPMs >  ", cpm.filt, " ..."))
    
    # 2025-05-04 Previously we were doing the following line isexpr = Matrix::rowMeans( cpm(counts, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)) >= cpm.filt # These are Aquino2023 settings: Check Methods > Pseudobulk estimation, normalization and batch correction
    # 1. Create DGEList
    dge <- edgeR::DGEList(counts = counts)
    # 2. Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    # 3. Calculate CPM with normalization
    cpm_matrix <- edgeR::cpm(counts, normalized.lib.sizes = TRUE) 
    # 4. Log2 transform and filter
    log2_cpm <- log2(cpm_matrix + 1)
    keep_genes <- rowMeans(cpm_matrix) >= cpm.filt
    print(paste0("N features expressed: ", table(keep_genes)["TRUE"], " / ", length(keep_genes), " ..."))
    
    # Output list
    list("raw_counts" = counts [ keep_genes, ], 
         "cpm_counts" = cpm_matrix [ keep_genes, ], 
         "log2_cpm_counts" = log2_cpm [ keep_genes, ]
    )
  }
}

cpm.norm <- function(counts, log.scale = FALSE){
  
  if(class(counts) == "list"){lapply(counts, function(x) cpm.norm( x, log.scale ) )
  }else{
  print(paste0("CPM normalization log=", log.scale))
  #/ make the DGEList:
  y <- DGEList(counts)
  
  #/ calculate TMM normalization factors:
  y <- calcNormFactors(y)
  
  #/ get the normalized counts:
  cpms <- cpm(y, log=log.scale)
  
  return(as.data.frame(cpms))
  }
}

##### RUN ##### 

## PARAMS
counts <- opt$counts
cpm.filt <- opt$cpm.filt
out_file <- opt$out_file

# 0. Read input file 
counts <- readRDS(counts)

# 1. Compute CPM normalization & Filter
# counts <- filter_by_cpm_counts( counts, cpm.filt = cpm.filt, frac.sample.filt = 0.4)
counts.l <- filter_by_mean_cpm_counts( counts, cpm.filt = cpm.filt ) # 2025-04-11 Switching to filtering by mean 

# 3. Save
## 3.0 Raw counts
saveRDS(counts.l[["raw_counts"]], out_file)

## 3.1 CPM counts
out_file <- gsub("raw-counts.rds", "cpm-counts.rds", out_file)
saveRDS(counts.l[["cpm_counts"]], out_file)

## 2.3 Log2 CPM counts
out_file <- gsub("cpm-counts.rds", "log2-cpm-counts.rds", out_file)
saveRDS(counts.l[["log2_cpm_counts"]], out_file)

print("Files saved ...")
