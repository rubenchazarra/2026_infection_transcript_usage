#!/usr/bin/env Rscript

# This script generates boxplot visualizations for infection-associated DTU genes
# across immune cell types and viral stimulations (influenza A virus and SARS-CoV-2).
# In this particular case, for the ORMDL1 gene represented in Figure 2D of the paper
#
# For each gene, four plots are generated:
# i)   Boxplot of differences in equivalence class (EC) ratio abundances between conditions (infected − baseline)
# ii)  Boxplot of raw EC counts by Condition (e.g COV and Baseline)
# iii) Boxplot of relative EC counts by condition
# iv)  Stacked barplot of relative EC counts per donor

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
    default = "ORMDL1",
    type = 'character',
    help = 'Gene name of the gene to plot ' 
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

## 0 . Packages
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))

## 0. Functions
current_dir <- getwd()
setwd(current_dir)
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
# source("../../../00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R") # 2025-07-11 Loading this reassigns names to old functions
# source("../../../00_Scripts/05_DTU-viz-02_EC-coord-Coverage-plot.R")
setwd(current_dir)

## DEFINE THEME------------------------------------------
plot_guide <-  guides( fill = FALSE, 
                       alpha=guide_legend( override.aes = list(label = "", size = 1, ncol=1, fill = "black")))
## PROCESSING FUNCTIONS ------------------------------------------
ec.palette.general <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# 2025-11-13 New palette for Paper
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

other.ECs.col <- "#6c757d"

order.ec.df <- function(ec.df){
  # print("Ordering Equivalence classes by N of transcripts they represent ...")
  # Add N of tx per EC
  ec.df[["ec_length"]] <- unlist( lapply(ec.df[["ec"]], function(ec){ length(unlist(strsplit(as.character(ec), ","))) }) )
  # Order
  ec.df <- ec.df[ order(ec.df[["ec_length"]], decreasing = F), ]
  # Order numerically tx of length 1
  ec.df[ec.df$ec_length == 1, ] <- ec.df[ order(ec.df[ec.df$ec_length == 1, "transcript_name" ], decreasing = F), ]
  return(ec.df)
}

# edit.bold.for.DGE.DTE <- function(annot.df) { # 2025-11-13 This function is not editing DTE transcripts in bold
#   # if DGE, gene_name subtitle in bold
#   gene_name <- unique(annot.df$gene_name)
#   sub_tit <- if(any(annot.df$is.DGE)){ bquote(bold(.(gene_name))) }else{ gene_name }
#   # if DTE, transcript_name in bold 
#   tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE)
#   tr.vec.DTE <- c()
#   for(i in 1:length(tr.vec)){
#     label <- if( annot.df[i, "is.DTE"] ){ bquote(bold(.( tr.vec[i] ))) }else{ tr.vec[i] }
#     tr.vec.DTE <- c(tr.vec.DTE, label)
#   }
#   names(tr.vec.DTE) <- annot.df$transcript_name
#   # if other.ECs in annotation, color differently
#   ec.palette <- ec.palette.general[1 : length(tr.vec)]; if("other.ECs"%in% tr.vec){ ec.palette[match("other.ECs", tr.vec)] <- other.ECs.col}
#   return(list("sub_tit" = sub_tit, "tr.vec.DTE" = tr.vec.DTE, "ec.palette" = ec.palette))
# }

sum.excluded.ECs <- function(exp.i, ec.annot.df, meta_data, annot.cols.excluded){
  # Given an annotation lowly abundant ECs (that should not appear in the plots), aggreate their RAW, CPM or Relative abundance COUNTS
  exp.ex <- as.data.frame( t( colSums(exp.i[ec.annot.df$exclude_plot, intersect(colnames(exp.i), meta_data$Donor.ID_Condition) ], na.rm = T) )) # sum excluded ECs
  exp.ex <- cbind(annot.cols.excluded, exp.ex) # annotate collapsed ECs excluded
  exp.i <-  exp.i[!ec.annot.df$exclude_plot, ] # Not excluded ECs
  exp.i <- rbind(exp.i, exp.ex) # merge both
  return(exp.i)
}

sum.excluded.ECs.wrapper <- function(raw_exp.i, rel_exp.i, diff_exp.i, ec.annot.df, meta_data ){
  # Given an annotation lowly abundant ECs (that should not appear in the plots), aggreate their RAW, CPM or Relative abundance COUNTS across all counts
  # NOte "<<-" symbol below updates the values in the global environment without the need to return # "<<-" symbol was not being applied when running from a different script which sources the present script
  # 2025-07-14 I don't understand the previous note, where is the "<<" is it needed ? 
  annot.cols <- c("gene_name", "gene_id", "transcript_name", "ec")
  annot.cols.excluded <- rel_exp.i[1, annot.cols] ; annot.cols.excluded$transcript_name <- "other.ECs" ; annot.cols.excluded$ec <- "other.ECs" ; rownames(annot.cols.excluded) <- "other.ECs"
  ### 2.2.1 Raw counts 
  raw_exp.i <- sum.excluded.ECs(raw_exp.i, ec.annot.df, meta_data, annot.cols.excluded)
  ### 2.2.2 Rel abundances
  rel_exp.i <- sum.excluded.ECs(rel_exp.i, ec.annot.df, meta_data, annot.cols.excluded)
  ### 2.2.3 Diff counts 
  diff_exp.ex <- diff_exp.i[ ec.annot.df$exclude_plot, ]
  diff_exp.ex <- as.data.frame( t( colMeans(diff_exp.i[ec.annot.df$exclude_plot, intersect(colnames(diff_exp.i), meta_data[["Donor.ID"]]) ], na.rm = T) )) # collapse excluded
  diff_exp.ex <- cbind(annot.cols.excluded, diff_exp.ex) # annotate collapsed ECs excluded
  diff_exp.i <- diff_exp.i [ !ec.annot.df$exclude_plot, ] # not excluded ECs
  diff_exp.i <- rbind(diff_exp.i, diff_exp.ex)  # merge both
  
  ### 2.2.5 Edit EC annotation
  ec.annot.df.ex  <- ec.annot.df[ ec.annot.df$exclude_plot, ]
  ec.annot.df.ex <- cbind(annot.cols.excluded, data.frame("ec_length" = paste(ec.annot.df.ex$ec_length, collapse = ","), 
                                                          "is.DGE" = FALSE, #ec.annot.df$is.DGE[1],
                                                          "is.DTE" = FALSE,
                                                          "EC.mean.abundance" = sum(ec.annot.df.ex$EC.mean.abundance), 
                                                          "exclude_plot" = F))
  ec.annot.df <- ec.annot.df[ ! ec.annot.df$exclude_plot, ]
  ec.annot.df <- rbind(ec.annot.df, ec.annot.df.ex)
  
  out.l <- setNames(list(raw_exp.i, rel_exp.i, diff_exp.i, ec.annot.df), nm = c("raw_exp.i", "rel_exp.i", "diff_exp.i", "ec.annot.df"))
  return(out.l)
}

## INFECTION DTU PLOT FUNCTIONS ------------------------------------------
raw_exp.plt.fun  <- function(exp.i, tit, annot.df, y_lab){
  # Generate a boxplot where x =EC, y = Raw counts, and alpha = Condition
  ## Shorten transcript names, instead of ABC-1,ABC-2 --> ABC-1,2
  tr.vec = annot.df$transcript_name ; tr.vec = paste0(paste0(gene_name, "-"), gsub(paste0(gene_name, "-"), "", tr.vec)) ; names(tr.vec) = annot.df$transcript_name
  
  # vals = edit.bold.for.DGE.DTE(annot.df)
  # tr.vec = vals[["tr.vec.DTE"]]; ec.palette = vals[["ec.palette"]]
  exp.i$transcript_name <- factor( exp.i$transcript_name, levels = names(tr.vec))
  
  ggplot(exp.i, aes(x= transcript_name, y = count, alpha = Condition )) + 
    geom_boxplot( aes(fill = transcript_name), outliers = F, linewidth = 0.25) + 
    scale_alpha_discrete("Condition", range=c(1, 0.5)) +
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    scale_x_discrete(labels = tr.vec) + 
    labs(title = tit) + ylab(y_lab) + 
    plot_theme + theme(plot.title = element_blank(),  #element_text(size = 6), 
                       axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
    plot_guide 
}

rel_exp.plt.per_Donor.ID <- function(rel_exp, tit, annot.df, y_lab){
  # Generate a barplot PER DONOR faceted per ancestry with rel tx values
  ## Shorten transcript names, instead of ABC-1,ABC-2 --> ABC-1,2
  tr.vec = annot.df$transcript_name ; tr.vec = paste0(paste0(gene_name, "-"), gsub(paste0(gene_name, "-"), "", tr.vec)) ; names(tr.vec) = annot.df$transcript_name
  rel_exp$transcript_name <- factor( rel_exp$transcript_name, levels = names(tr.vec))
  
  ggplot(rel_exp, aes(x= Donor.ID, y = count, fill= transcript_name)) + 
    ylim(0, 1) + 
    geom_bar( position="stack", stat="identity") +  
    facet_grid(~ Condition, scales = "free", space = "free" ) + 
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    # scale_x_discrete(labels = tr.vec) + 
    labs(title = tit) + ylab(y_lab) + 
    plot_theme + 
    theme(plot.title = element_blank(),  #element_text(size = 6), 
          axis.title.x= element_blank(), axis.line.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
}

infectionDTU.diff.Condition.plt <- function(diff_exp, tit, annot.df, y_lab){
  # Plot diference rel abundance (flu  - NI)
  ## Shorten transcript names, instead of ABC-1,ABC-2 --> ABC-1,2
  tr.vec = annot.df$transcript_name ; tr.vec = paste0(paste0(gene_name, "-"), gsub(paste0(gene_name, "-"), "", tr.vec)) ; names(tr.vec) = annot.df$transcript_name
  diff_exp$transcript_name <- factor( diff_exp$transcript_name, levels = names(tr.vec))
  
  ggplot(diff_exp, aes(x= transcript_name, y = count)) + 
    geom_boxplot( aes(fill = transcript_name), outlier.colour = "grey50", outlier.size = 0.5, linewidth = 0.25) +
    geom_hline(yintercept = 0, color = "black", linetype="dotdash") + 
    # scale_fill_brewer(leg.tit, labels = tr.vec, palette = "Set3")  +  
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  +  
    scale_x_discrete(labels = tr.vec) + 
    labs(title = tit) + ylab(y_lab) + 
    plot_theme + theme(plot.title = element_blank(),  #element_text(size = 6), 
                       axis.title.x = element_blank(), 
                       axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) + 
    plot_guide
}

infectionDTU.diff.raw.rel.barplot.wrapper <- function(dtu.l, diff.l, raw.l, rel.l, ct.condition, meta_data, gene_name){
  # This function computes: 1) TCC Raw counts Boxplot faceted by POP, 2) TCC Rel counts Boxplot faceted by POP, 3) TCC Rel Counts as Stacked Barplot
  # 0. Get DTU and Count data
  dtu.df = dtu.l[[ ct.condition ]] ; dtu.df[is.na(dtu.df)] <- F # 2025-02-19 BUG, DTU genes not appearing in the DGE or DTE shoudl appear  as F, not as NA, edit in DA script
  diff_exp = diff.l[[ ct.condition ]] # not list (its IAV - NS)
  raw_exp = raw.l[[ ct.condition ]] # no longer list, now grouped
  rel_exp = rel.l[[ ct.condition ]] # no longer list, now grouped
  
  # remove annotation from diff_exp and rel_exp (we will reannotate later)
  diff_exp <- diff_exp[, intersect(colnames(diff_exp), meta_data[["Donor.ID"]])]
  rel_exp <- rel_exp[, intersect(colnames(rel_exp), meta_data[["Donor.ID_Condition"]])]
  
  ct <- gsub(".NS|.COV|.IAV", "", ct.condition)
  condition.pairwise <- gsub(paste0(ct,"."), "", ct.condition)
  condition.pairwise.vec.short <- c("IAV.NS" = "IAV-Baseline",  "COV.NS" = "COV-Baseline", "IAV.COV" = "IAV-COV") 
  condition.pairwise <- plyr::revalue(condition.pairwise, condition.pairwise.vec.short)
  condition_1 <- strsplit(condition.pairwise, "-")[[1]][1] # this can be improved, we only use it in diff plot !
  condition_2 <- strsplit(condition.pairwise, "-")[[1]][2]
  
  # B) Iterate across genes
  dtu.df <- dtu.df[dtu.df$gene_name %in% gene_name, ] ## 2025-11-13 Added this line to make script work for a single genes
  lapply( dtu.df$gene_id, function(gene_id){
    tit <- ct.condition
    idx <- match(gene_id, dtu.df[["gene_id"]])
    # 1. Annotate EC ids
    ec.vec <- unlist(strsplit(dtu.df[idx , "feature_ids" ], ";"))
    if( ! dtu.df[idx, "feature_ids.excluded"] == "NA" ) {  ec.vec <- c( ec.vec, unlist(strsplit(dtu.df[idx , "feature_ids.excluded" ], ";"))) }
    ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 4 )
    gene_name <- unique(ec.annot.df$gene_name)
    ec.annot.df <- order.ec.df( ec.annot.df )
    # # add DGE info --> 2025-05-10 This is to highlight DTE, DGE, labels in the plot. Given DGE, DTU info is not present across all results, we are setting to false for now 
    # ec.annot.df$is.DGE <- dtu.df[idx, "is.DGE"]
    # # add DTE info
    # ec.dte <- unlist(strsplit(dtu.df[idx, "DTE.which"], ";"))
    # ec.annot.df$is.DTE <- F
    # if(any(ec.dte %in% ec.annot.df$ec)) { ec.annot.df[match(ec.dte, ec.annot.df$ec), "is.DTE"] <- T }
    ec.annot.df$is.DGE <- F
    ec.annot.df$is.DTE <- F
    gene_name <- unique(ec.annot.df$gene_name)
    ecs <- ec.annot.df$ec
    print(paste0(ct.condition, " - ", idx, " - ", gene_name, " - ", gene_id))
    
    # 2. Subset counts
    ## 2.1 Annotate and Subset Raw Counts
    raw_exp.i <- raw_exp[ecs, ] # subset
    raw_exp.i <- cbind(annotate.equivalence.classess(ecs, tx2gene, n_cores = 4), raw_exp.i)
    ## 2.2 Annotate and Subset Rel Counts
    rel_exp.i <- rel_exp[ecs, ] # subset
    rel_exp.i <- cbind(annotate.equivalence.classess(ecs, tx2gene, n_cores = 4), rel_exp.i)
    ## 2.3 Subset Diff (c_1 - c_2 ) counts
    diff_exp.i <- diff_exp[ ecs, ] # subset
    diff_exp.i <- cbind(annotate.equivalence.classess(ecs, tx2gene, n_cores = 4), diff_exp.i)
    
    # 3. Annotate ECs based on abundance (e.g label lowly abundant ECs 
    annot.cols <- c("gene_id", "gene_name", "transcript_name", "ec")
    ec.mean.abundances <- data.frame("EC.mean.abundance" = rowMeans(rel_exp.i[, intersect(meta_data$Donor.ID_Condition, colnames(rel_exp.i)) ], na.rm = T))
    ec.annot.df <- cbind(ec.annot.df, ec.mean.abundances [ ec.annot.df$ec, , drop = F ] )
    if( nrow(ec.annot.df) > 6 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.abundance < 0.05
    } else { ec.annot.df$exclude_plot <- F }
    ## 3.2 If lowly abundant EC present, collapse into "other.ECs" label
    if( length( ec.annot.df$exclude_plot[ec.annot.df$exclude_plot] ) > 1 ){ 
      out.l <- sum.excluded.ECs.wrapper(raw_exp.i, rel_exp.i, diff_exp.i, ec.annot.df, meta_data )
      raw_exp.i <- out.l[["raw_exp.i"]] ; rel_exp.i <- out.l[["rel_exp.i"]] ; diff_exp.i <- out.l[["diff_exp.i"]] ; ec.annot.df <- out.l[["ec.annot.df"]]
    }
    
    # 4. Plots
    tit <- paste0(tit, " - ", gene_name)
    ## 4.1 Raw expression plot 
    raw_exp.i <- reshape2::melt( raw_exp.i, value.name = "count", variable.name = "Donor.ID_Condition") # melt
    raw_exp.i$Condition <- gsub("..*_", "", raw_exp.i$Donor.ID_Condition) # add Condition
    raw_exp.i$Condition <- plyr::revalue(raw_exp.i$Condition, condition.vec.short)
    raw_exp.i$Condition <- factor(raw_exp.i$Condition, levels = c(condition_1, condition_2)) # put Baseline first in the plot
    raw_exp.plt <- raw_exp.plt.fun(raw_exp.i, tit = tit, ec.annot.df, y_lab = "Raw Counts")
    
    ### 4.2.1 BOXPLOT 
    rel_exp.i <- reshape2::melt( rel_exp.i, value.name = "count", variable.name = "Donor.ID_Condition") # melt
    rel_exp.i$Condition <- gsub("..*_", "", rel_exp.i$Donor.ID_Condition) # add Condition
    rel_exp.i$Condition <- plyr::revalue(rel_exp.i$Condition, condition.vec.short)
    rel_exp.i$Condition <- factor(rel_exp.i$Condition, levels = c(condition_1, condition_2)) # put Baseline first in the plot
    rel_exp.plt.1 <- raw_exp.plt.fun(rel_exp.i, tit = tit, ec.annot.df, y_lab = "EC Ratios")
    
    ### 4.2.2 Stacked BARPLOT per Donor.ID
    rel_exp.i$Donor.ID <- gsub("_.*", "", rel_exp.i$Donor.ID_Condition) # add Donor.ID
    rel_exp.plt.2 <- rel_exp.plt.per_Donor.ID(rel_exp.i, tit, ec.annot.df, y_lab = "EC Ratios")
    
    ## 4.3. Difference (c_1 - c_2 ) Plot
    diff_exp.i <- reshape2::melt( diff_exp.i, value.name = "count", variable.name = "Donor.ID" )
    y_lab <- paste0("EC Ratios ( ", condition_1, " - ", condition_2, " ) ")
    diff_exp.plt <- infectionDTU.diff.Condition.plt(diff_exp.i, tit, ec.annot.df, y_lab = y_lab)
    
    return(list(diff_exp.plt, raw_exp.plt, rel_exp.plt.1, rel_exp.plt.2))
  })
}

#### PARAMS
ct_condition <- opt$ct_condition
gene_name <- opt$gene_name
out_file <- opt$out_file

## RUN ------------------------------------------

# 1. Load data
## 1.1 popDTU data
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/01_infectionDTU/02_DA/01_Processed-Data/05_Filt-1-5-0.8-0.01/01_infectionDTU-filtered-list.rds"
dtu.l <- readRDS(file.name)
### 1.1.1 Select lineage 
dtu.l <- dtu.l$lineage
dtu.l <- lapply(dtu.l, function(ct) ct[c("COV.NS", "IAV.NS")])
## 1.2 Metadata
meta_data <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
## 1.3 Tx2gene
tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", header = T, sep = "\t")
## 1.4 Raw Counts
raw.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/04_pseudobulk/05_split-Condition/01_ECs/01_lineage/", pattern  = "group-Condition.rds")
## 1.5 Rel Counts
rel.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/02_infectionDTU-1-5-0.8-0.01/01_lineage/", pattern  = "group-Condition-rel_abs.rds")
## 1.6 Diff rel Counts
diff.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/02_infectionDTU-1-5-0.8-0.01/01_lineage/", pattern  = "group-Condition-rel_abs-diff-Condition.rds")
print("Data read ...")

# 2 PROCESS
## 2.1 Group RAW and REL counts by Condition 
# raw.l <- group_by_Condition_Aquino( raw.l )
# rel.l <- group_by_Condition_Aquino( rel.l )

## 2.2 Subset RAW Counts to transcripts from genes that are DTU ( for speed up )
# raw.l <- subset.rel.df(dtu.l, raw.l)
### 2.1.2 Subset Rel Counts
# rel.l <- subset.rel.df(dtu.l, rel.l)

## 2.2 Unlist 1 level
dtu.l <- unlist( dtu.l, recursive = F)
diff.l <- unlist( diff.l, recursive = F)
raw.l <- unlist( raw.l, recursive = F)
rel.l <- unlist( rel.l, recursive = F)

## 2.3 Select lineage_condition to plot
dtu.l <- dtu.l[ ct_condition ]
diff.l <- diff.l[ ct_condition ]
raw.l <- raw.l[ ct_condition ]
rel.l <- rel.l[ ct_condition ]

## 2.4 Select EC switch
# dtu.l <- lapply(dtu.l, function(x) x[x$is.EC.switch == "TRUE", ])

## 2.5 Fitler empty 
dtu.l <- dtu.l[ unlist(lapply(dtu.l, function(x) nrow(x) > 0)) ] 

## 2.6 Additional selection (e.g plot only EC switches)
dtu.l <- lapply(dtu.l, function(x) x[x$is.transcript.switch, ])

# 3. Plots
# A) Iterate across cell types.Condition
plt.l <- lapply( names(dtu.l), function( ct.condition ) infectionDTU.diff.raw.rel.barplot.wrapper(dtu.l = dtu.l, diff.l = diff.l, raw.l = raw.l, rel.l = rel.l, ct.condition, meta_data, gene_name) )
names(plt.l) <- names(dtu.l)
# fix names
# names(plt.l) <- gsub(".IAV", "_IAV", names(plt.l))
# names(plt.l) <- gsub(".COV", "_COV", names(plt.l))                     


## 4. Save plots
# path <- "03_infectionDTU-Viz/05_Filt-1-5-0.8-0.01/01_all/"
inch_to_cm <- 2.54
# lapply(names(plt.l), function(ct.condition){
#   print(paste0("Saving plots for ...", ct.condition))
#   file.name <- paste0(path, "01_infectionDTU-all-", ct.condition, "-Viz-01_BoxPlot.pdf")
pdf(out_file, height = 4/inch_to_cm, width = 6/inch_to_cm)
plt.l
dev.off()
print("Plots saved ...")
# })
