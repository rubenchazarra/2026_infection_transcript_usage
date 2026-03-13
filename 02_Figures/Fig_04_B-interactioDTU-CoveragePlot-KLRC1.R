#!/usr/bin/env Rscript

# FUNCTIONS
current_dir <- getwd()
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-02_EC-coord-Coverage-plot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
setwd(current_dir)
ec.palette.general <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# 2025-11-13 New palette for Paper
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# 1. Load data
## 1.1 popDTU data
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/03_interactionDTU/02_DA/01_Processed-Data/06_celltype-fractions/01_interactionDTU-filtered-list.rds"
dtu.l <- readRDS(file.name)
dtu.l = dtu.l$lineage
## 1.2 Tx2gene
tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", header = T, sep = "\t")
## 1.3 GTF 
gtf <- as.data.frame(rtracklayer::import("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_gencode.v42.annotation-filt-by-bulkRNA-seq.gtf"))
## 1.4 All ECs (Equivalence classes)
all.ecs <- readRDS("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/03_Merge-split-cts/01_tcc/01_lineage/00_equivalence-classes-lineages.rds")

# 2 PROCESS
## 2.2 Unlist 1 level
dtu.l <- unlist(dtu.l, recursive = F)
## 2.3 Select first N DTU hits
dtu.l <- select.first.N(dtu.l, sort.col = "fdr", N = Inf) # select all 
## 2.4 Remove empty 
dtu.l <- dtu.l[unlist(lapply(dtu.l, nrow)) > 0 ]

## 2.5 Select individual gene in individual lineage-condition 
ct.condition <- "T.CD8.IAV.NS"
gene_name <- "KLRC1"
dtu.l <- dtu.l[ct.condition]
dtu.l <- lapply(dtu.l, function(x) x[x$gene_name %in% gene_name, ])

# RUN
track.l <- lapply( names(dtu.l), function(ct.condition) interactionDTU.Tx.Cov.EqClass.track.wrapper(dtu.l, ct.condition ))
names(track.l) <- names(dtu.l)

# SAVE
suppressPackageStartupMessages(require(Gviz))
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/03_interactionDTU/02_DA/03_interactionDTU-Viz/06_celltype-fractions/02_CoveragePlots/02_indiv/01_Plots/T.CD8_IAV.NS_KLRC1_02_CoveragePlot-02_colors.pdf"  

pdf(file.name, width = 6.5, height = 6.5)
for (ct_name in names(track.l)) { # Iterate across cell types
  print(paste0(" ========== ", ct_name, " ========== "))
  for( gene_name in names(track.l [[ct_name]] )){ # Iterate across genes
    print(gene_name)
    track <- track.l[[ct_name]][[gene_name]]
    plotTracks(unlist(track, recursive = F),
               thinBoxFeature = c("utr"), # determines which feature should appear thin
               transcriptAnnotation = "symbol", 
               main = paste(ct_name, " - ", gene_name), 
               cex.main = 0.75, 
               cex = 1.5, 
               legend = T)
  }
}
dev.off()