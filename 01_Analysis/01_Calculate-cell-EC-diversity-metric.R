#!/usr/bin/env Rscript 

# 2025-07-14
# This script takes a Seurat object containing equivalence class (EC) counts
# for each immune lineage (e.g. CD4+ T, CD8+ T, monocytes, NK, and B cells),
# including cells from baseline, influenza A virus (IAV), and SARS-CoV-2 (COV)
# conditions, and performs the following steps:
#
# 1) Subsets the EC count matrix to retain the top 75% of cells with the highest
#    UMI counts, and downsamples each cell to a common minimum UMI count
#    (i.e. all retained cells are normalized to the same number of UMIs).
#
# 2) Transforms the EC × cell count matrix into a gene × cell matrix, where each
#    entry contains the number of ECs expressed (count > 0) for a given gene in
#    a given cell, and saves this object (output 1.1).
#
# 3) Calculates “cell EC diversity”, defined as the mean number of expressed ECs
#    (count > 0) per gene in each cell, restricting the analysis to genes with
#    at least two annotated transcripts (output 1.2).

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--seu"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Seurat object containing EC x cell counts'
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output Seurat object with N of exp ECs per gene per cell (matrix size: gene x cell)' 
  ), 
  make_option(
    c("-c", "--cell_metric"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output cell metric (matrix size: gene x cell)' 
  )#, 
  # make_option( # finally not included in publication
  #   c("-g", "--gene_metric"),
  #   action = "store",
  #   default = NA,
  #   type = 'character',
  #   help = 'Output gene metric' 
  # )
)

opt <- parse_args(OptionParser(option_list=option_list))

### PACKAGES ### 
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(scuttle))

### FUNCTIONS ### 
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R")

create.Seurat.with.n.ecs.exp.per.gene <- function(seu, count_filter, cell_n_filter){
  ## Input Seurat object of ECs x cells, filter ECs, calculate N of ECs expresed per gene per cell, and output Seurat object with genes x cells with N of Ecs per gene per cell expressed
  if(class(seu) == "list"){ lapply(seu, function(seu) create.Seurat.with.n.ecs.exp.per.gene(seu, count_filter, cell_n_filter))
  }else{
      print(dim(seu))
      meta <- seu@meta.data
      assays <- names(seu@assays)
      assay.l <- lapply(assays, function(x){ # iterate across each of the assays present in Seurat object
        mat = seu@assays[[x]]$counts
        # Filter ECs from matrix 
        keep.ecs = Matrix::rowSums(mat > count_filter) > cell_n_filter
        print(paste0("Count filter: ", count_filter, " | N of cells filter: ", cell_n_filter))
        print(table(keep.ecs))
        mat <- mat[ keep.ecs, ]
        # Annotate
        ec.annot = annotate.equivalence.classess(ec.vec = rownames(mat), tx2gene =tx2gene, n_cores = 112)
        mat@x[] <- 1
        gene.mat <-  SparseArray::rowsum( mat, group = ec.annot$gene_id )
        print("Summed ...")
        gene.mat <- gene.mat[ intersect(rownames(gene.mat), genes.more.1.tx), ]
        gene.mat
      })
      names(assay.l) <- assays
      # Create Seurat with RNA assay
      seu.gene <- Seurat::CreateSeuratObject(counts = assay.l[[assays[1]]], meta.data = meta, assay = assays[1] )
      
      # Add SCT assay
      if(length(assays) > 1){
        seu.gene[[assays[2]]] <- Seurat::CreateAssayObject(counts = assay.l[[assays[1]]])
      }
      
      print("Seurat generated")
      seu.gene
    }
}

calculate.EC.diversity.per.cell <- function(seu){
  # Compute Mean N of EC expressed per gene per cell across genes (one value per cell)
  if(class(seu) == "list"){ lapply(seu, calculate.EC.diversity.per.cell)
  }else{
    assays <- names(seu@assays)
    assay.l <- lapply(assays, function(x){
      mat = seu@assays[[x]]$counts
      print(dim(mat))
      # In mat we already have stored the N of ECs expressed (count > 0) per gene per cell
      n_ecs <- colSums(mat)
      n_genes <- colSums(mat > 0 ) 
      mean.n.ec.gene.cell <- n_ecs / n_genes 
      df <- seu@meta.data
      df[["mean_n_ec_gene_cell"]] <- mean.n.ec.gene.cell
      df[["n_ecs"]] <- n_ecs
      df[["n_genes"]] <- n_genes
      df
    })
    names(assay.l) <- assays
    assay.l
  }
}

# calculate.EC.diversity.per.gene <- function(seu){ # finally not used in publication
#   # Compute Mean N of EC expressed per gene per cell across cells (one value per gene)
#   ## Calculate for all cells and also splitting per population, this may be of interest later
#   if(class(seu) == "list"){ lapply(seu, calculate.EC.diversity.per.gene)
#   }else{
#     mat = seu@assays$RNA$counts
#     print(dim(mat))
#     df <- data.frame(gene_id = rownames(mat), 
#                      "mean_n_ec.all_cells" = rowMeans(mat),
#                      "mean_n_ec.EUB" = rowMeans(mat[, seu$POP %in% "EUB"]) ,
#                      "mean_n_ec.AFB" = rowMeans(mat[, seu$POP %in% "AFB"])
#     )
#     df.exp <- data.frame("mean_n_ec.all_cells.expressed" = apply(mat, 1, function(x) mean(x[x > 0])),
#                          "mean_n_ec.EUB.expressed" = apply(mat[, seu$POP %in% "EUB"], 1, function(x) mean(x[x > 0])),
#                          "mean_n_ec.AFB.expressed" = apply(mat[, seu$POP %in% "AFB"], 1, function(x) mean(x[x > 0]))
#     )
#     cbind(df, df.exp)
#   }
# }

# PARAMS
seu <- opt$seu
in_file <- seu
out_file <- opt$out_file
cell_metric <- opt$cell_metric
# gene_metric <- opt$gene_metric

# 1. Read Seurat
seu <- readRDS(seu)
print("File read ...")

## 1.2 Consider genes with > 1 Tx annotate, which will have > 1 EC
tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", sep = "\t", header = T)
tx.count <- tx2gene %>% group_by(gene_id, gene_name) %>% summarise (n_transcript_id = n_distinct(transcript_id), n_transcript_name = n_distinct(transcript_name))
genes.more.1.tx <- tx.count[ tx.count$n_transcript_id > 1, "gene_id", drop = T]

# 2. Process

## 2.1 NO downsampling
### 2.1.1 Create Seurat with N of ECs expressed /gene x cell matrix
seu.l <- Seurat::SplitObject(seu, "Condition") # if we don't split, the function below fails for large cell types with error Error in .Call2("C_rowsum_dgCMatrix", x, group, length(ugroup), na.rm,  : too many groups (matrix of sums will be too big)"
seu.gene <- create.Seurat.with.n.ecs.exp.per.gene(seu.l, count_filter = 1, cell_n_filter = 1) # We set cell_n_filter=1 to exclude ECs expressed in a single-cell
saveRDS(seu.gene, out_file)
### 2.1.1 Calculate Cell metric of EC diversity (Mean N of EC per gene per cell across genes, one value per cell)
cell.metric <- calculate.EC.diversity.per.cell(seu.gene)
saveRDS(cell.metric, cell_metric)
# ## 2.1.2 Calculate Gene metric of EC diversity (Mean N of EC per gene per cell across cells, one value per gene)
# gene.metric <- calculate.EC.diversity.per.gene(seu.gene)
# saveRDS(gene.metric, gene_metric)
print("EC diversity cell metric calculated with NO downsampling ...")

## 2.2 Downsampling to equal N of UMI counts per cell
## Rationale --> To do this, we downsample single-cell UMI count matrix to a common minimum of UMIs, avoid the common minimum UMI count being too low, we exclude the 25% of cells with less UMIs
### 2.1.0 Subsample seurat object removing 25% of cells with less UMIs
n_umis_per_cell <- colSums(seu@assays$RNA$counts) 
keep_cells <- n_umis_per_cell >= quantile(n_umis_per_cell, 0.25) # & n_count <= quantile(n_count, 0.75) # 2025-08-20: We retain cells with 75% top UMI counts
seu <- seu[, keep_cells]
### 2.1.1 Downsample matrix to equal N of UMIs per cell
mat <- seu@assays$RNA$counts
mat <- mat[rowSums(mat) > 0, ] # Remove zero count features after cell filterings
n_umis_per_cell <- colSums(mat)
p <- min(n_umis_per_cell) / n_umis_per_cell
mat.d <- scuttle::downsampleMatrix(mat, prop = p, bycol = TRUE)
n_umis_downsampled <- colSums(mat.d)
print(summary(n_umis_downsampled)) # Check N UMI distribution after downsampling
### 2.1.2 Create downsampled Seurat object
seu.d <- Seurat::CreateSeuratObject(counts = mat.d, meta.data = seu@meta.data )
seu.d@meta.data[["N_UMIS_downsampled"]] <- n_umis_downsampled
### 2.1.3 Save downsampled Seurat object 
out_file.d <- gsub(".rds", "-downsample-equal-nCount-per-cell.rds", out_file) 
saveRDS(seu.d, out_file.d)

### 2.1.1 Create Seurat with N of ECs expressed /gene x cell matrix 
seu.d.l <- Seurat::SplitObject(seu.d, "Condition") # if we don't split, the function below fails for large cell types with error Error in .Call2("C_rowsum_dgCMatrix", x, group, length(ugroup), na.rm,  : too many groups (matrix of sums will be too big)"
seu.gene.d <- create.Seurat.with.n.ecs.exp.per.gene(seu.d.l, count_filter = 1, cell_n_filter = 1) # We set cell_n_filter=1 to exclude ECs expressed in a single-cell
### 2.1.4 Calculate Cell metric of EC diversity downsampled (Mean N of EC per gene per cell across genes, one value per cell)
cell.metric <- calculate.EC.diversity.per.cell(seu.gene.d)
cell_metric <- gsub(".rds", "-downsample-equal-nCount-per-cell.rds", cell_metric) 
saveRDS(cell.metric, cell_metric)

# ## 2.1.5 Calculate Gene metric of EC diversity downsampled (Mean N of EC per gene per cell across cells, one value per gene)
# gene.metric <- calculate.EC.diversity.per.gene(seu.gene.d)
# gene_metric <- gsub(".rds", "-downsample-equal-nCount-per-cell.rds", gene_metric) 
# saveRDS(gene.metric, gene_metric)
print("EC diversity cell metric calculated with WITH downsampling ...")

print("All file saved ...")