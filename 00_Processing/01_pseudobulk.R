#!/usr/bin/env RScript

# Calculate pseudo-bulk by aggregating counts per Sample.ID for each immune lineage (e.g. CD4+ T, CD8+ T, monocytes, NK, and B cells)

suppressPackageStartupMessages(require("optparse"))
option_list = list(
  make_option(
    c("-i", "--seu"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Seurat object to Pseudo-bulk'
  ),
  make_option(
    c("-a", "--assay"),
    action = "store",
    default = "RNA",
   type = 'character',
    help = 'Assay to pseudobulk from: either RNA (unnnormalized) or SCT (normalized)'
  ),
  make_option(
    c("-f", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# FUNCTIONS 
pseudobulk <- function(seu){
  if(class(seu) == "list"){lapply(seu, pseudobulk)
    }else{
      # Create sample ID. Note each sample in Aquino et al.(2023) has replicates, hence we need to create a now covariate identifying cells for each sample. 
      seu@meta.data[["Sample.ID"]] <- paste0( seu@meta.data[["Library"]], "_", seu@meta.data[["POP"]],  "_", seu@meta.data[["Donor.ID"]], "_", seu@meta.data[["Condition"]] )
      # Get count matrix 
      mat  <- GetAssayData( seu, assay = assay, slot = "counts" )
      # Calculate pseudobulk by summing counts
      meta_data <- seu@meta.data
      ids <- meta_data[["Sample.ID"]]
      unique_ID_list <- as.list(unique(meta_data[["Sample.ID"]]))
      pseudobulk <- as.data.frame( pblapply( unique_ID_list, FUN = function(x){ sparse_Sums ( mat [, ids == x, drop = FALSE], rowSums = TRUE ) } ) )
      # add names
      colnames(pseudobulk) <- unique_ID_list
      rownames(pseudobulk) <- rownames(mat)
      return(pseudobulk)
  }
}

# PACKAGES
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(textTinyR))
suppressPackageStartupMessages(require(pbapply))

# PARAMS 
seu <- opt$seu
assay <- opt$assay
out_file <- opt$out_file


# 1. Read RDS objects
seu <- readRDS(seu)
print("Seurat object read ...")

# 2. Pseudobulk
pb <- pseudobulk(seu)

# 3. Save iobject list as unique file
saveRDS(pb, file = out_file)
print("File saved ...")
