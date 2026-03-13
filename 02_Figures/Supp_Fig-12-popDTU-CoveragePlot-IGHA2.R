#!/usr/bin/env Rscript

# This script generates a locus-level visualization for genes identified in the
# population differential transcript usage (DTU) analysis.
#
# For each selected gene, the script:
# 1) retrieves the annotated transcript structures corresponding to the
#    equivalence classes (ECs) retained in the DTU model;
# 2) computes the genomic coordinates associated with each EC based on the
#    transcript regions represented by that EC;
# 3) displays the EC-specific genomic segments over the transcript models;
# 4) computes mean read coverage across genomic positions separately for
#    EUB (European) and AFB (AFrican) samples; and
# 5) plots the coverage tracks together with the transcript and EC annotations
#    to facilitate interpretation of the DTU signal.
#
# The resulting figure integrates:
# i)   transcript models for the gene of interest,
# ii)  genomic coordinates of the modeled ECs,
# iii) mean coverage in infected samples, and
# iv)  mean coverage in baseline samples.
#
# This visualization is used to interpret how changes in EC relative abundance
# map onto transcript structure and read support across the locus.
# Corresponds to Figure 2d of the paper.

# Parameter values Rscript ../../../00_Scripts/03_popDTU-Viz-02_CoveragePlot_02_indiv.R -g IGHA2  -l B -c NS -o 01_Plots/  --height 12     --width 12

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-g", "--gene_name"),
    action = "store",
    default = NA,
    type = 'character',
    help = ' gene_name of gene to be plotted ' # IGHA2
  ), 
  make_option(
    c("-l", "--celltype"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Celltype ' # B
  ), 
  make_option(
    c("-c", "--condition"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Condition. Must be one of:NS (Baseline), IAV (influenza A), COV (SARS-CoV2)' # NS
  ), 
  make_option(
    c("-o", "--out_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output_path ' 
  ), 
  make_option(
    c("-t", "--height"),
    action = "store",
    default = 5,
    type = 'numeric',
    help = 'out file plot height ' 
  ), 
  make_option(
    c("-w", "--width"),
    action = "store",
    default = 7,
    type = 'numeric',
    help = 'out file plot width ' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list, add_help_option=T))

# PACKAGES 
suppressPackageStartupMessages(require(Gviz))
# FUNCTIONS
current_dir <- getwd()
common.path <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/"
source(paste0(common.path, "01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R"))
source(paste0(common.path, "01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R"))
source(paste0(common.path, "02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R"))
source(paste0(common.path, "02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-02_EC-coord-Coverage-plot.R"))
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
setwd(current_dir)
# 2025-11-13 New palette for Paper
ec.palette.general <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# 0. Params
gene_name <- opt$gene_name
celltype <- opt$celltype
condition <- opt$condition
height <- opt$height
width <- opt$width
out_path <- opt$out_path
out_file <- paste0(out_path, "/",  celltype, "-", condition, "-", gene_name, "-popDTU-02_CoveragePlot.pdf")
inch.to.cm <- 2.54 

# 1. Load data
## 1.1 popDTU data
path <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/"
file.name <- paste0(path, "03_Differential/02_DTU/02_DTU/02_popDTU/02_DA/01_Processed-Data/06_celltype-fractions/01_popDTU-filtered-list.rds")
dtu.l <- readRDS(file.name)
dtu.l <- dtu.l$lineage
## 2025-08-07 I am keeping the list structure to avoid having duplicate functions single-gene plot generation
dtu.l <- unlist(dtu.l, recursive = F)
dtu.l <- dtu.l[ paste0(celltype, ".", condition ) ] # select cell type and condition 
dtu.l <- lapply(dtu.l, function(x) x[x$gene_name %in% gene_name, ])

## Tx2gene
tx2gene <- paste0(path, "/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv")
tx2gene <- read.table(tx2gene, header = T, sep = "\t")
## GTF 
gtf <-  paste0(path, "/01_Quantif/00_IndexData/02_gencode.v42.annotation-filt-by-bulkRNA-seq.gtf")
gtf <- as.data.frame(rtracklayer::import(gtf))
## set of all Equivalence classes # 205-0'4-16 I think this should be cell type specific
ecs <-  paste0(path,"02_Processing/03_Merge-split-cts/01_tcc/01_lineage/00_equivalence-classes-lineages.rds")
all.ecs <- readRDS(ecs)

# 3. RUN 
track.l <- lapply( names(dtu.l), function(ct.condition) popDTU.Tx.Cov.EqClass.track.wrapper(dtu.l = dtu.l, ct.condition ))
names(track.l) <- names(dtu.l)

# 4. SAVE
pdf(out_file, height= height / inch.to.cm, width = width / inch.to.cm)
for (ct_name in names(track.l)) { # Iterate across cell types
  print(paste0(" ========== ", ct_name, " ========== "))
  for( gene_name in names(track.l [[ct_name]] )){ # Iterate across genes
    print(gene_name)
    track <- track.l[[ct_name]][[gene_name]]
    plotTracks(unlist(track, recursive = F),
               thinBoxFeature = c("utr"), # determines which feature should appear thin
               transcriptAnnotation = "symbol", 
               main = paste(ct_name, " - ", gene_name), 
               cex.main = 0.1, 
               cex = 0.1, 
               legend = T)
  }
}
dev.off()