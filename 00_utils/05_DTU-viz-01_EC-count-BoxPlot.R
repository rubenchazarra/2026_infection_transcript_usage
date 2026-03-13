#!/usr/bin/env Rscript

# 2025-07-14 Now plots are in scripts of the individual comparisons: infectionDTU, popDTU, interaction DTU

# 2025-02-17 Functions for visualization of ECs as boxplots
## 1) Between infection status (infection DTU)
## 2) Between populations (population DTU)
## 3) Between populations upon infection (interaction DTU)
## 4) Across all cell types for each of the previous ones

suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(Gviz))

##  0. Common functions ------------------------------------------
group_by_Condition_pairwise <- function(counts.l){
  print("Restructuring list to pairwise Condition comparison ...")
  # Group list of Condition count dataframes ready for pairwise comparison: "COV-NS", "IAV-NS", "COV-IAV"
  groups <- list("IAV.NS" = c("IAV", "NS"), "COV.NS" = c("COV", "NS"), "IAV.COV" = c("IAV", "COV"))
  
  if(class(counts.l[[1]]) == "list"){ # if object inside list is list, iterate. This is for when we have data split by POP
    lapply(counts.l, group_by_Condition_pairwise)
  }else{
    counts.l <- lapply(groups, function(x) {
      l <- list(counts.l[[ x[1] ]], counts.l[[ x[2] ]])
      setNames(l, nm = x)
    })
    names(counts.l) <- names(groups)
    return(counts.l)
  }
}

ec.palette.general <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")
other.ECs.col <- "#6c757d"

edit.bold.for.DGE.DTE <- function(annot.df) {
  # if DGE, gene_name subtitle in bold
  gene_name <- unique(annot.df$gene_name)
  sub_tit <- if(any(annot.df$is.DGE)){ bquote(bold(.(gene_name))) }else{ gene_name }
  # if DTE, transcript_name in bold 
  tr.vec = gsub("((?:[^,]*,){1}[^,]*),", "\\1\n", annot.df$transcript_name, perl = TRUE)
  tr.vec.DTE <- c()
  for(i in 1:length(tr.vec)){
    label <- if( annot.df[i, "is.DTE"] ){ bquote(bold(.( tr.vec[i] ))) }else{ tr.vec[i] }
    tr.vec.DTE <- c(tr.vec.DTE, label)
  }
  names(tr.vec.DTE) <- annot.df$transcript_name
  # if other.ECs in annotation, color differently
  ec.palette <- ec.palette.general[1 : length(tr.vec)]; if("other.ECs"%in% tr.vec){ ec.palette[match("other.ECs", tr.vec)] <- other.ECs.col}
  return(list("sub_tit" = sub_tit, "tr.vec.DTE" = tr.vec.DTE, "ec.palette" = ec.palette))
}

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

##  0. Common plots:  RAW EC Counts, Relative abundance EC ratios # Valid for both: infection DTU & population DTU ------------------------------------------
##  1. infection DTU DIFF, RAW, REL EXP PLOTS (TCC/TX)  Valid for both: infection DTU & population DTU------------------------------------------

popDTU.sum.excluded.ECs <- function(exp, ec.annot.df, annot.cols){
  # Given an annotation lowly abundant ECs (that should not appear in the plots), aggreate their RAW, CPM or Relative abundance COUNTS
  # sum excluded
  annot.df <- exp [1 , annot.cols ]
  annot.df$ec <- "other.ECs" ; annot.df$transcript_name <- "other.ECs"
  exp.excluded <- as.data.frame( t( colSums(exp [ ec.annot.df$exclude_plot, ! colnames(exp) %in% annot.cols], na.rm = T) )) # collapse excluded
  exp.excluded <- cbind( annot.df, exp.excluded )
  # not excluded
  exp <- exp[ ! ec.annot.df$exclude_plot , ]
  # merge both
  exp <-  rbind( exp, exp.excluded ) 
  return(exp)
}

###  1.2. Interaction Plots ------------------------------------------ # 2025-07-14 I think now in interaction plot
# diff.condition.int.plt.fun <- function(exp.i, tit, sub_tit, annot.df, leg.tit = "TCC", y_lab = "Diff rel abundance (cond_1 − cond_2)"){ # Same function as diff.flu_ni.plt.fun (above) only different in facet_wrap(~POP)
#   # Plot diference rel abundance (flu  - NI)
#   tr.vec = gsub(",", "\n", annot.df$transcript_name) ; names(tr.vec) = annot.df$transcript_name
#   exp.i[["transcript_name"]] <- factor(exp.i[["transcript_name"]], levels = annot.df$transcript_name)
#   ec.palette <- ec.palette[1 : length(tr.vec)]; if("other.ECs"%in% tr.vec){ ec.palette[match("other.ECs", tr.vec)] <- other.ECs.col}
#   
#   ggplot(exp.i, aes(x= transcript_name, y = count)) + 
#     geom_boxplot( aes(fill = transcript_name), outlier.colour = "grey50", outlier.size = 0.5) +
#     geom_hline(yintercept = 0, color = "red", linetype="dotdash") + 
#     facet_wrap( ~ POP ) + 
#     # scale_fill_brewer(leg.tit, labels = tr.vec, palette = "Set3")  +
#     scale_fill_manual(leg.tit, labels = tr.vec, values = ec.palette)  + 
#     scale_x_discrete( labels = tr.vec ) + 
#     labs(title = tit, subtitle = sub_tit) + ylab(y_lab) + 
#     theme(panel.background = element_blank(), plot.background = element_blank(), 
#           axis.title.x = element_blank(), axis.text.x = element_blank(),
#           axis.line.x = element_line(),
#           #strip.text.x = element_text(size = 6, angle = 90), 
#           axis.line.y = element_line(), 
#           legend.text = element_text(size = 6), 
#           legend.spacing.y = unit(0.25, 'cm')) + 
#     guides(fill = guide_legend(byrow = TRUE))
# }
# 
# raw_exp.int.plt.fun  <- function(exp.i, tit, sub_tit, annot.df, leg.tit = "TCC", grid_col = "infection_status"){
#   # Generate a boxplot per tx faceted per ancestry with raw tx values
#   tr.vec = gsub(",", "\n", annot.df$transcript_name) ; names(tr.vec) = annot.df$transcript_name
#   exp.i[["transcript_name"]] <- factor(exp.i[["transcript_name"]], levels = annot.df$transcript_name)
#   ec.palette <- ec.palette[1 : length(tr.vec)]; if("other.ECs"%in% tr.vec){ ec.palette[match("other.ECs", tr.vec)] <- other.ECs.col}
#   
#   ggplot(exp.i, aes(x= transcript_name, y = count )) + 
#     geom_boxplot( aes(fill = transcript_name), outlier.colour = "grey50", outlier.size = 0.5) + 
#     facet_wrap(~ get(grid_col), scal) + 
#     # scale_fill_brewer(leg.tit, labels = tr.vec, palette = "Set3") +
#     scale_fill_manual(leg.tit, labels = tr.vec, values = ec.palette)  + 
#     scale_x_discrete(labels = tr.vec) + 
#     labs(title = tit, subtitle = sub_tit) + ylab("CPM Counts") + 
#     theme(panel.background = element_blank(), plot.background = element_blank(), 
#           axis.title.x = element_blank(), axis.text.x = element_blank(), 
#           axis.line.x = element_line(),
#           #strip.text.x = element_text(size = 6, angle = 90), 
#           axis.line.y = element_line(), 
#           legend.text = element_text(size = 6), 
#           legend.spacing.y = unit(0.25, 'cm')) +  
#     guides(fill = guide_legend(byrow = TRUE))
# }
# 
