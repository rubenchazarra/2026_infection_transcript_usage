#!/usr/bin/env Rscript

# This script generates a locus-level visualization for genes identified in the
# differential transcript usage (DTU) analysis.
#
# For each selected gene, the script:
# 1) retrieves the annotated transcript structures corresponding to the
#    equivalence classes (ECs) retained in the DTU model;
# 2) computes the genomic coordinates associated with each EC based on the
#    transcript regions represented by that EC;
# 3) displays the EC-specific genomic segments over the transcript models;
# 4) computes mean read coverage across genomic positions separately for
#    infected and baseline samples; and
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

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--ct_condition"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Lineage_Condition combination to plot DTU cases from. Must follow the syntax e.g "B.IAV" (with dot in the middle) ' 
  ), 
  make_option(
    c("-g", "--gene_name"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Gene_name to generate the plot from' 
  ), 
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output PDF file name' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Gviz))

# FUNCTIONS
current_dir <- getwd()
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-02_EC-coord-Coverage-plot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
setwd(current_dir)
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# PARAMS 
ct_condition <- opt$ct_condition
gene_name <- opt$gene_name
out_file <- opt$out_file

# 1. Load data
## 1.1 infectionDTU data --> Load list for visualization
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/01_infectionDTU/02_DA/01_Processed-Data/05_Filt-1-5-0.8-0.01/01_infectionDTU-filtered-list.rds" 
dtu.l <- readRDS(file.name)
dtu.l = dtu.l$lineage
dtu.l <- lapply(dtu.l, function(ct) ct[c("IAV.NS", "COV.NS")])
## 1.2 Rel Counts --> 2025-07-29 What is this for ? 
rel.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/02_infectionDTU-1-5-0.8-0.01/01_lineage/", pattern  = "group-Condition-rel_abs.rds")
## 1.3 Tx2gene
tx2gene <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv"
tx2gene <- read.table(tx2gene, header = T, sep = "\t")
## 1.4 GTF 
gtf <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_gencode.v42.annotation-filt-by-bulkRNA-seq.gtf"
gtf <- as.data.frame(rtracklayer::import(gtf))
## 1.5 ECs -- set of all Equivalence classes
ecs <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/03_Merge-split-cts/01_tcc/01_lineage/00_equivalence-classes-lineages.rds" # All ECs including library specific
all.ecs <- readRDS(ecs)
## Intersection of ECs --> 2025-07-29 When using the ECs that are common to all libraries, in many cases the gene tested for DTU has no ECs common to all libraries, hence cannot plot it !
# ecs <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/03_Merge-split-cts/01_tcc/01_lineage/PBMCs-ECs-common-to-all-libs-EUB-AFB.tsv" # Common ECs to all libraries
# all.ecs <- read.table(ecs, sep = "\t", header = T)
## 1.6 Metadata
meta_data <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
print("Data reads ...")
# 2 PROCESS

## 2.2 Unlist 1 level 
dtu.l <- unlist(dtu.l, recursive = F)
rel.l <- unlist(rel.l, recursive = F)

## 2.4 Subset to isoform switches
# print("Subsetting to Transcript Switch ...")
# dtu.l <- lapply(dtu.l, function(x) x[x$is.transcript.switch, ])

## 2.1 Select first N DTU hits
dtu.l <- select.first.N(dtu.l, sort.col = "fdr", N = Inf)

## 2.3 Select lineage_condition to plot
dtu.l <- dtu.l[ ct_condition ]
rel.l <- rel.l[ ct_condition ]

## 2.4 Subset to gene_name to plot
dtu.l[[ ct_condition ]] <- dtu.l[[ct_condition]][ match(gene_name, dtu.l[[ct_condition]]$gene_name), ]
rel.l[[ ct_condition ]] <- rel.l[[ct_condition]][ rel.l[[ct_condition]]$gene_name %in% gene_name, ]

# 3. RUN
track.l <- lapply( names(dtu.l), function(ct.condition) infectionDTU.Tx.Cov.EqClass.track.wrapper(dtu.l, rel.l, ct.condition ))
names(track.l) <- names(dtu.l)

# 4. SAVE
inch_to_cm <- 2.54
  
pdf(out_file, width=15/inch_to_cm, height = 12/inch_to_cm)
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
               # cex = 1.5, 
               legend = T)
  }
}
dev.off()