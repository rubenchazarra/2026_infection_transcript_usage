#!/usr/bin/env Rscript 

# This script contains functions called by differentdifferential gene expression scripts. 
# In the present study we run DEG across different covariates: infection DEG, population DEG and infection*population interaction DEG.

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(variancePartition)) # Dream tutorial in: https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/variancePartition/inst/doc/dream.html
suppressPackageStartupMessages(require(BiocParallel))
suppressPackageStartupMessages(require(stats))

# Common functions 
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_Scripts/01_DEA_common_functions.R")

load.Aquino.metadata <- function(){
  print("Loading Aquino2023 sample metadata ...")
  # sample_meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv"
  # 2025-05-04 Updated metadata, now one value per donor, and "Library" covariate represents a combination of the libraries of the 2 replicate samples of each Donor.ID_Condition combination
  sample_meta_data <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_per_Donor.ID_cell_fractions.tsv"
  sample_meta_data <- read.table(sample_meta_data, header = TRUE, sep = "\t")
  sample_meta_data$POP <- factor(sample_meta_data$POP, levels = c("EUB", "AFB"))
  sample_meta_data$Condition <- factor(sample_meta_data$Condition, levels = c("NS", "COV", "IAV"))
  return(sample_meta_data)
}

dge.limma.dream.wrapper  <- function(counts, sample_col, sample_meta_data, covariates, param ){
  # Iterativelly run DGE upon counts 
  if(class(counts) == "list"){ lapply(counts, function(counts) dge.limma.dream.wrapper(counts, sample_col, sample_meta_data, covariates, param) )
  }else{
    print("Running limma dream DGE...")
    # 1. Edit METADATA ( include case samples )
    print(dim(counts))
    sample_cols <- intersect( colnames(counts), sample_meta_data[[sample_col]] )
    meta_data <- sample_meta_data[ match(sample_cols, sample_meta_data[[sample_col]]), ]
    rownames(meta_data) <- meta_data[[sample_col]]
    # 2. Edit model.
    ## 2.1 Sanity check: if all values in one metadata column are the same, remove covariate
    # covariates <- covariates[ sapply(meta_data[, covariates], function(x) length(unique(x)) > 1) ]
    ## 2.2 Add random factors
    if("Run" %in% covariates) { print("Adding 'Run' as random effect ..."); covariates[ grep("Run", covariates) ] <- "(1|Run)" }
    if("Library" %in% covariates) { print("Adding 'Library' as random effect ..."); covariates[ grep("Library", covariates) ] <- "(1|Library )" }
    if("Donor.ID" %in% covariates) { print("Adding 'Donor.ID' as random effect ..."); covariates[ grep("Donor.ID", covariates) ] <- "(1|Donor.ID)" }
    ## 2.3 Model 
    model <- as.formula(paste(" ~  ", paste(c(covariates), collapse = " + ")))
    print( model )
    # 1. Normalize by library size
    dge <- DGEList(counts = counts)
    ## 2. Normalize by library size
    dge <- calcNormFactors(dge) # I think this has already been performed
    # 3. Estimate weights using linear mixed model of dream
    vobjDream <- voomWithDreamWeights( counts = dge, model, meta_data, span = 0.01, BPPARAM = param)
    # 4. Estimate weights using linear mixed model of dream
    vobjDream <- dream( vobjDream, model, meta_data, BPPARAM = param)
    # 5. run eBayes
    vobjDream = variancePartition::eBayes(vobjDream)
    # 6. Collect output
    vout <- collect_vfit_output(vobjDream)
    ## add to list
    return(vout)
  }
}


