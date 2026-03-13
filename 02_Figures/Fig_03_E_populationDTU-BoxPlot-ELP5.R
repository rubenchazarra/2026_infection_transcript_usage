#!/usr/bin/env Rscript

# This script generates boxplot visualizations for population -associated DTU genes
# across immune cell types and stimulations conditions (Baseline, influenza A, SARS-CoV-2).
# In this particular case, for the IGHA2 gene represented in Supplementary Figure 12B of the paper
#
# For each gene, three plots are generated:
# i)  Boxplot of raw EC counts by Condition (e.g COV and Baseline)
# ii) Boxplot of relative EC counts by condition
# iii)  Stacked barplot of relative EC counts per donor (Supplementary Figure 12.B)

# Parameter values are Rscript  ../../../00_Scripts/03_popDTU-Viz-01_BP_02_indiv.R -g ELP5 -l NK -c NS -o 01_Plots/ --width 5 --height 4.5


suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-g", "--gene_name"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'gene_name of gene to be plotted ' # ELP5
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
    help = 'Condition ' # NS
  ), 
  make_option(
    c("-o", "--out_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output_path ' 
  ), 
  make_option(
    c("-h", "--height"),
    action = "store",
    default = 5,
    type = 'numeric',
    help = 'out file plot height '  # 5
  ), 
  make_option(
    c("-w", "--width"),
    action = "store",
    default = 7,
    type = 'numeric',
    help = 'out file plot width ' # 4.5
  )
)

opt <- parse_args(OptionParser(option_list=option_list,  add_help_option=FALSE))

# FUNCTIONS
current_dir <- getwd()
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-02_EC-coord-Coverage-plot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
setwd(current_dir)

# 2025-11-13 New palette for Paper
ec.palette.general <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

plot.pre.processing <- function(exp, ec.annot.df){
  # Data processing for plot, including adding POP to Donor.ID_Conditions
  ## Melt dataframe
  exp.i <- reshape2::melt( exp, variable.name = "Donor.ID_Condition", value.name = "count")
  ## Merge metadata
  exp.i <- merge(exp.i, unique(meta_data[, c("Donor.ID_Condition", "POP")]), by ="Donor.ID_Condition")
  return(exp.i)
}

raw_exp.plt.fun  <- function(exp.i, tit, sub_tit, annot.df, y_lab ){
  # Generate a boxplot per tx faceted per ancestry with raw tx values
  ## Shorten transcript names, instead of ABC-1,ABC-2 --> ABC-1,2
  tr.vec = annot.df$transcript_name ; tr.vec = paste0(paste0(gene_name, "-"), gsub(paste0(gene_name, "-"), "", tr.vec)) ; names(tr.vec) = annot.df$transcript_name
  # tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE); names(tr.vec) = annot.df$transcript_name
  exp.i$transcript_name <- factor(exp.i$transcript_name, levels = names(tr.vec))
  
  ggplot(exp.i, aes(x= transcript_name, y = count, alpha = POP)) + 
    geom_boxplot( aes(fill = transcript_name), outliers = F, linewidth=0.15) + # currently ignoring outliers, to include them add this line # outlier.colour = "grey50", outlier.size = 0.5) + 
    scale_alpha_discrete("Population", range=c(0.3, 1)) +
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    scale_x_discrete(labels = tr.vec) + 
    labs(title = tit, subtitle = sub_tit) + ylab(y_lab) + 
    plot_theme + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=0, hjust = 0.5, vjust=0.5)) + 
    guides( fill = FALSE, alpha=guide_legend( override.aes = list(label = "", size = 1, ncol=1, fill = "grey20")))
}

order.Donor.ID.by.EC.rel <- function(rel_exp){
  # Order Donor.ID_Condition by ascending order of rel exp of most highly expressed tx
  max.rel_exp <- rel_exp %>% group_by(transcript_name) %>% summarise(mean_rel_exp = mean(count, na.rm = T))
  max.tx <- max.rel_exp[which.max(max.rel_exp$mean_rel_exp), "transcript_name", drop = T]
  # order by rel exp of most expressed tx
  max.rel_exp <- rel_exp[rel_exp$transcript_name %in% max.tx, ]
  max.rel_exp <- max.rel_exp[order(max.rel_exp$count), ]
  # factor Donor.ID_Condition
  rel_exp$Donor.ID_Condition <- factor(rel_exp$Donor.ID_Condition,
                                       levels = max.rel_exp$Donor.ID_Condition)
  return(rel_exp)
}

rel_exp.plt.per_Donor.ID <- function(rel_exp.i, tit, sub_tit, annot.df, y_lab){
  # Generate a stacked barplot per donor faceted per ancestry with rel tx values
  ## Shorten transcript names, instead of ABC-1,ABC-2 --> ABC-1,2
  tr.vec = annot.df$transcript_name ; tr.vec = paste0(paste0(gene_name, "-"), gsub(paste0(gene_name, "-"), "", tr.vec)) ; names(tr.vec) = annot.df$transcript_name
  # tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE); names(tr.vec) = annot.df$transcript_name
  rel_exp.i$transcript_name <- factor(rel_exp.i$transcript_name, levels = names(tr.vec))
  
  ggplot(rel_exp.i, aes(x= Donor.ID_Condition, y = count, fill= transcript_name)) + 
    geom_bar( position="stack", stat="identity") + 
    facet_grid(~ POP, scales = "free_x" ) + 
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    labs(title = tit, subtitle = sub_tit) + ylab(y_lab) + 
    plot_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    guides(fill = guide_legend(byrow = TRUE))
}

# PARAMS 
gene_name <- opt$gene_name
celltype <- opt$celltype
condition <- opt$condition
celltype.condition <- paste0(celltype, ".", condition)
height <- opt$height
width <- opt$width
out_path <- opt$out_path
out_file <- paste0(out_path, "/",  celltype, "-", condition, "-", gene_name, "-popDTU-01_BoxPlot.pdf")

# RUN 
print(paste0("Plots for gene=  ", gene_name, "  | lineage=  ", celltype, "  | Condition=  ", condition, " "))

# 1. Load data
## 1.1 popDTU data
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/02_popDTU/02_DA/01_Processed-Data/06_celltype-fractions/01_popDTU-filtered-list.rds"
dtu.l <- readRDS(file.name)
## 1.2 Metadata
meta_data <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
## 1.3 Tx2gene
tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", header = T, sep = "\t")
## 1.4 Raw Counts ( split_Condition )
raw.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/04_pseudobulk/05_split-Condition/01_ECs/01_lineage", pattern  = "raw-counts-avg-reps-split-Condition.rds")
## 1.5 Rel Counts ( split_Condition )
rel.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/01_popDTU-1-5-0.8-0.01/01_lineage", pattern  = "split-Condition-rel_abs.rds")
print("Files read ...")

# 2. Process
## 2.1 Subset counts to select input: gene, celltype, condition
dtu.l <- dtu.l$lineage
dtu <- dtu.l[[celltype]][[condition]]
dtu <- dtu[dtu$gene_name %in% gene_name, ]
raw <- raw.l[[celltype]][[condition]]
ec.annot.df <- annotate.equivalence.classess(ec.vec = rownames(raw), tx2gene, n_cores = 10)
ec.annot.df <- ec.annot.df[ec.annot.df$gene_name %in% gene_name, ]
ecs <- ec.annot.df$ec
raw <- raw[ecs, ]
# annotate raw
raw <- merge( ec.annot.df, raw, by.x = "ec", by.y = 0 )
rel <- rel.l[[celltype]][[condition]]
rel <- rel[ rel$gene_name %in% gene_name, ] 

## 2.2. Annotate EC ids
idx <- match(gene_name, dtu[["gene_name"]])
ec.vec <- unlist(strsplit(dtu[idx ,"feature_ids" ], ";"))
ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 4 )
# gene_name <- unique(ec.annot.df$gene_name)  #not needed

## 2.3 Annotate ECs based on abundance (e.g label lowly abundant ECs)
annot_cols <- c("gene_id", "gene_name", "transcript_name", "ec")
sample_cols <- intersect(colnames(rel), meta_data[["Donor.ID_Condition"]] )
rel <- rel[, c(annot_cols, sample_cols) ]
ec.mean.abundances <- data.frame("EC.mean.abundance" = rowMeans(rel[, !colnames(rel) %in% annot_cols], na.rm = T))
ec.annot.df$EC.mean.abundance <- ec.mean.abundances[["EC.mean.abundance"]] [ match(ec.annot.df$ec, rownames(ec.mean.abundances)) ]
rel.thres <- 0 # 0.05 # Do not exclude any EC for now
ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.abundance < rel.thres
# Order ECs (ECs representing 1 transcript first, later ECs representing >1 transcript)
ec.annot.df <- order.ec.df(ec.annot.df)
# 2025-03-21 Management of cases where lowly abundant ECs are collapsed to "other.ECs" in standby. TODO --> ensure same colnames in rel and raw for later plot pre.processing, and ensure correct ofder of other.ECs
# if( nrow(ec.annot.df) > 5 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.abundance < 0.05 } else { ec.annot.df$exclude_plot <- F }
# ## 3.2 If lowly abundant EC present, collapse into "other.ECs" label
# if( any(ec.annot.df$exclude_plot) ){ 
#   rel <- popDTU.sum.excluded.ECs(rel, ec.annot.df, annot_cols )
#   raw <- popDTU.sum.excluded.ECs(raw, ec.annot.df, annot_cols = c("gene_id", "gene_name", "transcript_name", "ec"))
#   # must re-order ec.annot.df
#   ec.annot.df <- raw[, c("gene_id", "gene_name", "transcript_name", "ec")]
# }

# 3. Plots
## 3.1 Raw expression
tit <- paste0(celltype, "_", condition) # No title will be in the plot
raw <- plot.pre.processing(raw, ec.annot.df)
raw_exp.plt <- raw_exp.plt.fun(raw, tit, sub_tit = gene_name, ec.annot.df, y_lab = "Raw Counts")

## 3.2 Rel expression Plots
rel <- plot.pre.processing(rel, ec.annot.df)
### A) BOXPLOT 
rel_exp.plt.1 <- raw_exp.plt.fun(rel, tit, sub_tit = gene_name, ec.annot.df, y_lab = "EC ratios")
### B) Stacked Boxplot
rel <- order.Donor.ID.by.EC.rel( rel )
rel_exp.plt.2 <- rel_exp.plt.per_Donor.ID(rel, tit, sub_tit  = gene_name, ec.annot.df, y_lab = "EC ratios")

plot.l <- list(raw_exp.plt, rel_exp.plt.1, rel_exp.plt.2)

# 3.3 Save plots
inch_to_cm <- 2.54
pdf(out_file, width = width/inch_to_cm, height = height/inch_to_cm )
plot.l
dev.off()

# out_file <- gsub(".pdf", ".png", out_file)
# png(out_file, width = width/inch_to_cm, height = height/inch_to_cm, bg = "white")
# plot.l
# dev.off()