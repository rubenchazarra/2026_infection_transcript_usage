## 0. Packages
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(dplyr))

## 0. Functions
current_dir <- getwd()
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/02_DTU/00_Scripts/00_MANTA-functions.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/00_Scripts/05_DTU-viz-01_EC-count-BoxPlot.R")
source("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/03_Aesthetics-repo.R")
setwd(current_dir)

ec.palette.general <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

# 2025-11-13 New palette for Paper
ec.palette.general <- c("#355070", "#b56576", "#6d597a", "#eaac8b", "#e56b6f", # Immune lineage colors
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")

condition.color.vec.short <- c("Baseline"="#969696", "IAV"= "black", "COV" = "black") 
line_width <- 0.25

subset.rel.df <- function(dtu.df, rel.df){
  if(class(dtu.df) == "list" & class(rel.df) == "list"){ mapply(subset.rel.df, dtu.df, rel.df, SIMPLIFY = F) # iterate if list
  }else{
    ecs <- unlist(strsplit(dtu.df[["feature_ids"]], ";"))
    ecs <- c( ecs, unlist(strsplit( dtu.df[["feature_ids.excluded"]], ";")))
    rel.df [ ecs , ]
  }
}

add.meta_data <- function(exp, meta_col = "Donor.ID_Condition"){
  # Data processing for plot, including adding POP to Donor.ID_Conditions
  get.from.meta <- c("POP", "Condition")
  if(!"Donor.ID" %in% colnames(exp)) get.from.meta <- c( get.from.meta, "POP" )
  meta.i <- meta_data[ match(exp[[meta_col]], meta_data[[ meta_col ]]), c("POP", "Condition") ]
  exp <- cbind(exp, meta.i)
  # Factor Condition
  exp$Condition <- plyr::revalue(exp$Condition, condition.vec.short)
  exp$Condition <- factor(exp$Condition, levels = c(condition_2, condition_1)) # put Baseline first in the plot
  exp$POP <- factor(exp$POP, levels = c("AFB", "EUB"))
  # Create POP_Condition
  exp$POP_Condition <- paste0(exp$POP, "_", exp$Condition )
  exp$POP_Condition <- factor(exp$POP_Condition, levels = c("AFB_Baseline", "AFB_IAV", "AFB_COV", "EUB_Baseline", "EUB_IAV", "EUB_COV"))
  return(exp)
}


raw_exp.plt.fun.interactionDTU_02  <- function(exp.i, tit, annot.df, y_lab ){
  # Generate a boxplot per tx faceted per ancestry with raw tx values
  tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE); names(tr.vec) = annot.df$transcript_name
  exp.i$transcript_name <- factor(exp.i$transcript_name, levels = names(tr.vec))
  if(length(tr.vec) > 8){plot_theme <- plot_theme +  theme(axis.text.x = element_text(size = 4), legend.title= element_text(size=5), legend.text = element_text(size=4), legend.key.size = unit(0.25, "cm")) }
  
  ggplot(exp.i, aes(x= transcript_name, y = count )) + 
    geom_boxplot( aes(fill = transcript_name, alpha = POP, color = Condition), outliers = T, outlier.size = 0.25, linewidth = 0.25) + # currently ignoring outliers, to include them add this line # outlier.colour = "grey50", outlier.size = 0.5) + 
    scale_alpha_discrete("Population", range=c(0.3, 1)) +
    # facet_grid(~ Condition) + 
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    scale_color_manual("Condition", values = condition.color.vec.short) + 
    scale_x_discrete(labels = tr.vec) + 
    labs(title = tit) + ylab(y_lab) + 
    plot_theme + theme(plot.title = element_blank(), 
                       axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
    guides(alpha = guide_legend( override.aes = list(label = "", size = 1, ncol=1, fill = "grey")), fill = FALSE )
}


order.Donor.ID.by.EC.rel <- function(rel_exp){
  # Order Donor.ID_Condition by ascending order of rel exp of most highly expressed tx
  max.rel_exp <- rel_exp %>% group_by(transcript_name) %>% summarise(mean_rel_exp = mean(count, na.rm = T))
  max.tx <- max.rel_exp[which.max(max.rel_exp$mean_rel_exp), "transcript_name", drop = T]
  # order by rel exp of most expressed tx
  max.rel_exp <- rel_exp[rel_exp$transcript_name %in% max.tx, ]
  max.rel_exp <- max.rel_exp[order(max.rel_exp$count), ]
  # factor Donor.ID_Condition
  rel_exp$Donor.ID_Condition <- factor(rel_exp$Donor.ID_Condition, levels = max.rel_exp$Donor.ID_Condition)
  return(rel_exp)
}

rel_exp.plt.per_Donor.ID <- function(rel_exp.i, tit, annot.df, y_lab){
  # Generate a stacked barplot per donor faceted per ancestry with rel tx values
  # tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE); names(tr.vec) = annot.df$transcript_name
  tr.vec = gsub(",", "\n", annot.df$transcript_name, perl = F); names(tr.vec) = annot.df$transcript_name # Replace "," in trasncript name with "\n"
  rel_exp.i$transcript_name <- factor(rel_exp.i$transcript_name, levels = names(tr.vec))
  if(length(tr.vec) > 8){plot_theme <- plot_theme +  theme(axis.text.x = element_text(size = 4), legend.title= element_text(size=5), legend.text = element_text(size=4), legend.key.size = unit(0.25, "cm")) }
  
  ggplot(rel_exp.i, aes(x= Donor.ID_Condition, y = count, fill= transcript_name)) + 
    geom_bar( position="stack", stat="identity") + 
    facet_grid(~ POP_Condition, scales = "free_x" ) + 
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    labs(title = tit) + ylab(y_lab) + 
    plot_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    guides(fill = guide_legend(byrow = TRUE))
}

diff.condition.int.plt.fun <- function(exp.i, tit, annot.df,  y_lab ){ 
  # Plot diference rel abundance (flu  - NI) alphaing by POPulation
  # tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE); names(tr.vec) = annot.df$transcript_name
  tr.vec = gsub(",", "\n", annot.df$transcript_name, perl = F); names(tr.vec) = annot.df$transcript_name # Replace "," in trasncript name with "\n"
  exp.i$transcript_name <- factor(exp.i$transcript_name, levels = names(tr.vec))
  if(length(tr.vec) > 8){ plot_theme <- plot_theme +  theme(axis.text.x = element_text(size = 4), legend.text = element_text(size=4)) }
  
  ggplot(exp.i, aes(x= transcript_name, y = count, alpha = POP)) + 
    geom_hline(yintercept = 0, color = "black", linetype="dotdash") + 
    geom_boxplot( aes(fill = transcript_name), outlier.colour = "grey50", outlier.size = 0.25, linewidth = 0.25) +
    scale_fill_manual("EC", labels = tr.vec, values = ec.palette.general)  + 
    scale_x_discrete( labels = tr.vec ) +
    scale_alpha_discrete("Population", range=c(0.3, 1)) +
    labs(title = tit) + ylab(y_lab) + 
    plot_theme +  theme(plot.title = element_blank(), 
                        axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
    guides( fill = FALSE,  alpha=guide_legend( override.aes = list(label = "", size = 1, ncol=1, fill = "grey")))
}

barplot.wrapper.interactionDTU.top <- function(dtu.l, raw.l, rel.l, diff.l, ct.condition){
  # This function computes: 1) EC Raw counts Boxplot faceted by POP, 2) EC Rel counts Boxplot faceted by POP, 3) EC Rel Counts as Stacked Barplot
  
  # This function computes: 1) EC Raw counts Boxplot faceted by POP, 2) EC Rel counts Boxplot faceted by POP, 3) EC Rel Counts as Stacked Barplot
  # 0. Get DTU and Count data
  dtu.df = dtu.l[[ ct.condition ]] ; dtu.df[is.na(dtu.df)] <- F # 2025-02-19 BUG, DTU genes not appearing in the DGE or DTE shoudl appear  as F, not as NA, edit in DA script
  # Avoid plotting genes with too many ECs
  dtu.df <- dtu.df[dtu.df$n_features < 20, ]
  diff = diff.l[[ ct.condition ]] # not list (its IAV - NS)
  raw = raw.l[[ ct.condition ]] # list: IAV, NS
  rel = rel.l[[ ct.condition ]] # list: IAV, NS
  
  # Remove annotation from diff_exp and rel_exp (we will reannotate later)
  diff <- diff[, intersect(colnames(diff), meta_data[["Donor.ID"]])]
  rel <- rel [, intersect(colnames(rel), meta_data[["Donor.ID_Condition"]]) ]
  ct <- gsub(".IAV.NS|.COV.NS|.IAV.COV", "", ct.condition)
  condition.pairwise <- gsub(paste0(ct,"."), "", ct.condition)
  condition.pairwise.vec.short <- c("IAV.NS" = "IAV-Baseline",  "COV.NS" = "COV-Baseline", "IAV.COV" = "IAV-COV") 
  condition.pairwise <- plyr::revalue(condition.pairwise, condition.pairwise.vec.short)
  condition_1 <<- strsplit(condition.pairwise, "-")[[1]][1] # Set as global variables to be used in add.meta_data()
  condition_2 <<- strsplit(condition.pairwise, "-")[[1]][2]
  
  # B) Iterate across genes
  lapply( dtu.df$gene_id, function(gene_id){
    
    i <- match(gene_id, dtu.df[["gene_id"]])
    # 1. Annotate EC ids
    ec.vec <- unlist(strsplit(dtu.df[i , "feature_ids" ], ";"))
    if( ! is.na(dtu.df[i, "feature_ids.excluded"]) ) {  ec.vec <- c( ec.vec, unlist(strsplit(dtu.df[i , "feature_ids.excluded" ], ";"))) }
    ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 4 )
    gene_name <- unique(ec.annot.df$gene_name)
    ec.annot.df <- order.ec.df( ec.annot.df )
    # # add DGE info --> 2025-05-10 This is to highlight DTE, DGE, labels in the plot. Given DGE, DTU info is not present across all results, we are setting to false for now 
    # ec.annot.df$is.DGE <- dtu.df[i, "is.DGE"]
    # # add DTE info
    # ec.dte <- unlist(strsplit(dtu.df[i, "DTE.which"], ";"))
    # ec.annot.df$is.DTE <- F
    # if(any(ec.dte %in% ec.annot.df$ec)) { ec.annot.df[match(ec.dte, ec.annot.df$ec), "is.DTE"] <- T }
    ec.annot.df$is.DGE <- F
    ec.annot.df$is.DTE <- F
    gene_name <- unique(ec.annot.df$gene_name)
    ecs <- ec.annot.df$ec
    tit <- paste0(ct.condition, " - ", gene_name)
    print(paste0(ct.condition, " - ", i, " - ", gene_name, " - ", gene_id))
    
    # 2. Subset counts
    ## 2.1 Annotate and Subset Raw Counts
    raw.i <- raw [ ecs, ] # subset
    raw.i <-  cbind(annotate.equivalence.classess(rownames(raw.i), tx2gene, n_cores = 4), raw.i) # add annotation
    ## 2.1 Annotate and Subset Rel Counts
    rel.i <- rel [ ecs, ] # subset
    rel.i <-  cbind(annotate.equivalence.classess(rownames(rel.i), tx2gene, n_cores = 4), rel.i) # add annotation
    ## 2.3 Subset Diff (c_1 - c_2 ) counts
    diff.i <- diff[ ecs, ] # subset
    diff.i <- cbind(annotate.equivalence.classess(rownames(diff.i), tx2gene, n_cores = 4), diff.i)
    
    # 3. Annotate ECs based on abundance (e.g label lowly abundant ECs 
    # annot.cols <- c("gene_id", "gene_name", "transcript_name", "ec")
    # ec.mean.abundances <- lapply(rel_exp.i, function(x) data.frame("EC.mean.abundance" = rowMeans(x[, intersect(meta_data$Donor.ID_Condition, colnames(x)) ], na.rm = T)))
    # ec.mean.abundances <- setNames( cbind( ec.mean.abundances[[1]], ec.mean.abundances[[2]] ), nm = c("EC.mean.1", "EC.mean.2") ) # should match because we are working with the intersection of ECs in condition_1 and condition_2
    # ec.annot.df <- cbind(ec.annot.df, ec.mean.abundances [ ec.annot.df$ec, ] )
    # if( nrow(ec.annot.df) > 6 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.1 < 0.025 & ec.annot.df$EC.mean.2 < 0.025 
    # } else { ec.annot.df$exclude_plot <- F }
    
    ## 3.2 If lowly abundant EC present, collapse into "other.ECs" label
    # if( length( ec.annot.df$exclude_plot[ec.annot.df$exclude_plot] ) > 1 ){ 
    #   out.l <- sum.excluded.ECs.wrapper(raw_exp.i, rel_exp.i, diff_exp.i, ec.annot.df, meta_data )
    #   raw_exp.i <- out.l[["raw_exp.l"]] ; rel_exp.i <- out.l[["rel_exp.l"]] ; diff_exp.i <- out.l[["diff_exp"]] ; ec.annot.df <- out.l[["ec.annot.df"]]
    # }
    
    # 4. Plots
    ## 4.1 Raw expression plot
    raw.i <- reshape2::melt( raw.i, value.name = "count", variable.name = "Donor.ID_Condition") # melt
    raw.i <- add.meta_data( raw.i, meta_col = "Donor.ID_Condition")
    raw_exp.plt <- raw_exp.plt.fun.interactionDTU_02(raw.i, tit = tit, ec.annot.df, y_lab = "Raw Counts")
    
    ### 4.2.1 BOXPLOT 
    rel.i <- reshape2::melt( rel.i, value.name = "count", variable.name = "Donor.ID_Condition") # melt
    rel.i <- add.meta_data( rel.i, meta_col = "Donor.ID_Condition")
    rel_exp.plt.1 <- raw_exp.plt.fun.interactionDTU_02(rel.i, tit = tit, ec.annot.df, y_lab = "EC Ratios")
    
    ### 4.2.2 Stacked BARPLOT per Donor.ID
    rel.i = order.Donor.ID.by.EC.rel(rel.i)
    rel_exp.plt.2 <- rel_exp.plt.per_Donor.ID(rel.i, tit, ec.annot.df, y_lab = "EC Ratios")
    
    ## 4.3. Difference (c_1 - c_2 ) Plot
    diff.i <- reshape2::melt( diff.i, value.name = "count", variable.name = "Donor.ID" )
    diff.i <- add.meta_data( diff.i, meta_col = "Donor.ID")
    y_lab <- paste0("EC Ratios ( ", condition_1, " - ", condition_2, " ) ")
    diff_exp.plt <- diff.condition.int.plt.fun(diff.i, tit, ec.annot.df, y_lab = y_lab)
    
    return(list(diff_exp.plt, raw_exp.plt, rel_exp.plt.1, rel_exp.plt.2))
  })
}

#### RUN 
# 1. Load data
## 1.1 popDTU data
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/03_interactionDTU/02_DA/01_Processed-Data/06_celltype-fractions/01_interactionDTU-filtered-list.rds"
# file.name <- "01_Processed-Data/05_Filt-1-5-0.8-0.01/01_interactionDTU-filtered-list.rds"
dtu.l <- readRDS(file.name)
dtu.l <- dtu.l$lineage # select lineage
## 1.2 Metadata
meta_data <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
## 1.3 Tx2gene
tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", header = T, sep = "\t")
## 1.4 Raw Counts (group_Condition)
raw.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/02_Processing/04_pseudobulk/05_split-Condition/01_ECs/01_lineage/", pattern  = "group-Condition.rds")
## 1.5 Rel Counts (group_Condition)
rel.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/02_infectionDTU-1-5-0.8-0.01/01_lineage/", pattern  = "group-Condition-rel_abs.rds")
## 1.5 Rel Counts
diff.l <- read.RDS.files(path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/01_rel_abundances/02_infectionDTU-1-5-0.8-0.01/01_lineage/", pattern  = "diff-Condition.rds")
print("Files read ...")

# 2 PROCESS
### 2.1.1. Subset Raw Counts to transcripts from genes that are DTU ( for speed up )
raw.l <- subset.rel.df(dtu.l, raw.l)
### 2.1.2 Subset Rel Counts
rel.l <- subset.rel.df(dtu.l, rel.l)
### 2.1.2 Subset Rel Counts
diff.l <- subset.rel.df(dtu.l, diff.l)

## 2.2 Unlist 1 level
dtu.l <- unlist(dtu.l, recursive = F)
raw.l <- unlist(raw.l, recursive = F)
rel.l <- unlist(rel.l, recursive = F)
diff.l <- unlist(diff.l, recursive = F)

## Remove empty 
dtu.l <- dtu.l[unlist(lapply(dtu.l, nrow)) > 0 ]

## 2.3 Select first N DTU hits 
dtu.l <- select.first.N ( dtu.l, sort.col = "fdr", N = Inf ) # Selecting all

## 2.4 Select gene of interest 
ct.condition <- "T.CD8.IAV.NS"
gene_name <- "KLRC1"
dtu.l <- dtu.l[ct.condition]
dtu.l <- lapply(dtu.l, function(x) x[x$gene_name %in% gene_name, ])

# 3. Plots
# A) Iterate acrss cell types.Condition
plt.l <- lapply( names(dtu.l), function(ct.condition) barplot.wrapper.interactionDTU.top(dtu.l, raw.l, rel.l, diff.l, ct.condition)) 
#names(plt.l) <- names(dtu.l)[1]

## Save plots
file.name <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/02_DTU/02_DTU/03_interactionDTU/02_DA/03_interactionDTU-Viz/06_celltype-fractions/01_BoxPlots/02_indv/01_Plots/01_interactionDTU-01_BoxPlots-T.CD8.IAV.NS_KLRC1.pdf"
inch.to.cm <- 2.54
pdf(file.name, height = 4.5/inch.to.cm, width = 7.5/inch.to.cm)
plt.l
dev.off()
print("Plots saved ...")