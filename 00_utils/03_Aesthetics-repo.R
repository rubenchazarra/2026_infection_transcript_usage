 ## 2024-10-03 This script is to contain the Aesthetic directives of the Single-cell AS projec related to the Aquino2023 dataset

##  0. PACKAGES  ------------------------------------------
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(ComplexUpset))
##  0. LOAD RECURRENT DATA  ------------------------------------------
# tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/00_IndexData/02_GENCODEv42_filt-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", sep = "\t", header = T)


##  1. PLOT THEME  ------------------------------------------
plot_theme <- theme(panel.background = element_blank(), plot.background = element_blank(), # background
                    plot.title = element_blank(), plot.subtitle = element_blank(), # element_text(size=6),
                    axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), # axis title
                    axis.text.x = element_text(size = 5, angle = 0),
                    axis.text.y = element_text(size = 5, angle = 0),
                    axis.line.x = element_line(linewidth = 0.25), 
                    axis.line.y = element_line(linewidth = 0.25), # line
                    legend.title = element_text(size = 5, margin = margin(l=0.05, r=0.05, t=0.05, b=0.05, unit = "cm")),
                    legend.text = element_text(size = 5, margin = margin(l=0.05, r=0.05, t=0.05, b=0.05, unit = "cm")),
                    legend.key.size = unit(0.25, "cm"),
                    legend.key.height = unit(0.25, "cm"),
                    legend.key.spacing.x = unit(0.25, "cm"),
                    legend.key.width = unit(0.25, "cm"),
                    legend.spacing.y = unit(0.05, "cm"),
                    legend.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "cm"),
                    legend.box.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "cm"),
                    legend.box.spacing = unit(0, "cm"),
                    axis.ticks = element_line(linewidth = 0.25), # ticks
                    strip.background = element_rect(fill = NA, color = "white", size = 0.25), # facets
                    strip.placement = "outside", 
                    # strip.text = element_text(size = 7 ),
                    strip.text.x = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")), 
                    strip.text.y = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")), 
                    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), unit ="cm")
)

# 1. Cell count at the lineage (broad) and cell type (refined) level after all the preprocessing conducted
setwd(dirname(normalizePath(sys.frame(1)$ofile)))
# 1.2 Lineage color values
color.df <- data.frame(color = c( "T.CD4" = "#355070", "T.CD8" = "#b56576", "MONO" = "#6d597a" ,  "NK" = "#eaac8b", "B" =  "#e56b6f" ))
color.df$lineage = rownames(color.df)
lineage.color.vec <- setNames(color.df$color, nm = color.df$lineage)
lineage.df <- read.table( "01_metadata_exploration/03_Aquino2023-filtered-cell-count-01_lineage.txt", sep = "\t", header = T)
lineage.df <- merge(lineage.df, color.df)
lineage.df <- lineage.df[ order(lineage.df$n_cells, decreasing = T), ]
## 1.3 cell_type color
celltype.df <- read.table("01_metadata_exploration/03_Aquino2023-filtered-cell-count-02_celltype.txt", sep = "\t", header = T)
celltype.df <- merge(celltype.df, color.df, by = "lineage")
celltype.df <- celltype.df[ order(celltype.df$n_cells, decreasing = T), ]
celltype.color.vec <- setNames(celltype.df$color, nm = celltype.df$celltype)
celltype.lineage.vec <- setNames(celltype.df$lineage, nm = celltype.df$celltype)

# 2. Condition color values 
condition.vec <- c("NS" = "Baseline", "IAV" = "Influenza A", "COV" = "SARS-CoV2")
condition.vec.short <- c("NS" = "Baseline", "IAV" = "IAV", "COV" = "COV")
condition.color.vec <- c( "Baseline" = "#969696", "Influenza A" = "#75A9D1", "SARS-CoV2" = "#c1121f" )
condition.color.vec.short <- c( "Baseline" = "#969696", "IAV" = "#75A9D1", "COV" = "#c1121f" )

# 2.2. Pairwise Condition comparison color values 
condition.pairwise.vec <- c("IAV.NS" = "Influenza A - Baseline",  "COV.NS" = "SARS-CoV2 - Baseline", "IAV.COV" = "Influenza A - SARS-CoV2") 
condition.pairwise.vec.short <- c("IAV.NS" = "IAV-Baseline",  "COV.NS" = "COV-Baseline", "IAV.COV" = "IAV-COV") 

condition.pairwise.color.vec = c( "Influenza A - Baseline" = "#75A9D1", "SARS-CoV2 - Baseline" = "#c1121f", "Influenza A - SARS-CoV2" = "#c8b6ff")
condition.pairwise.color.vec.short = c( "IAV-Baseline" = "#75A9D1", "COV-Baseline" = "#c1121f",  "IAV-COV" = "#c8b6ff")


##  2. Plot N cells  ------------------------------------------

line_width <- 0.25
plot_theme.aes <- theme(panel.background = element_blank(), plot.background = element_blank(), # background
                    plot.title = element_blank(), plot.subtitle = element_blank(), # element_text(size=6),
                    # axis.title.x = element_blank(), axis.title.y = element_blank(), # axis title
                    axis.text.x = element_text(size = 5, angle = 0),
                    axis.text.y = element_text(size = 5, angle = 0),
                    axis.line.x = element_line(linewidth = line_width), 
                    axis.line.y = element_line(linewidth = line_width), # line
                    legend.title=element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")), 
                    legend.text = element_text(size = 5, margin = margin(r = 10, unit = "pt")),
                    legend.key.size = unit(0.25, "cm"), 
                    legend.margin=margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
                    legend.box.spacing = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
                    axis.ticks = element_line(linewidth = line_width), # ticks
                    strip.background = element_rect(fill = NA, color = "white", size = line_width), # facets
                    strip.placement = "outside", 
                    # strip.text = element_text(size = 7 ),
                    strip.text.x = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")), 
                    strip.text.y = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")), 
                    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
)

plot_guide <-  guides( fill = FALSE, # guide_legend(byrow = TRUE), 
                       alpha=guide_legend( override.aes = list(label = "", size = 1, ncol=1, fill = "grey")))

##  2.1 N of cells per Lineage  ------------------------------------------
suppressPackageStartupMessages(require(ggplot2))
# N cells
generate.n.cells.per.lineage <- function(){
  
  df <- read.table("01_metadata_exploration/03_Aquino2023-filtered-cell-count-01_lineage-Condition.txt")
  
  df$Condition <- plyr::revalue(df$Condition, condition.vec)
  df$Condition <- factor(df$Condition, levels = rev(condition.vec) )
  color.vec <- setNames(lineage.df$color, nm = lineage.df$lineage)
  df$lineage <- factor(df$lineage, levels = names(color.vec))
  
  # Plot 
  x_lab <- "N cells"
  y_lab <- "Cell type"
  
  df$label_color <- ifelse(df$lineage %in% c("T.CD4", "MONO"), "white", "black")
  # 
  n_cell.plt <- ggplot(data = df, aes(x = n_cells, y = Condition, fill = lineage)) +
    geom_bar(stat= "identity", linewidth=0.25, position = position_stack()) + 
    # Labels
    geom_text(data = df, aes(x = n_cells, label = n_cells),
              # color = ifelse(lineage %in% c("T.CD4", "MONO"), yes = "white", no = "black")), # this doesnt work evet!
              size = 3.5, hjust = 0.5, position = position_stack(vjust = .5)) +
    facet_grid(lineage ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100) * 1e3, labels = c("0", "25K", "50K", "75K", "100K")) + 
    scale_fill_manual( values = lineage.color.vec ) +
    xlab(x_lab) + ylab(y_lab) + 
    theme(plot.title = element_blank(), plot.subtitle = element_blank(),
          panel.background = element_blank(), plot.background = element_blank(),
          axis.line.x = element_line(), axis.line.y = element_line(),
          axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
          #  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.position = c(0.75, 0.2), legend.box = "vertical",
          legend.title=element_text(size=9), legend.text=element_text(size=8),
          strip.background = element_rect(fill = NA, color = "black", size = 1), 
          strip.placement = "outside") + 
    guides(fill = FALSE, color = FALSE ) 
  
  n_cell.plt
}

n_cell.lineage.plt <- generate.n.cells.per.lineage()


generate.n.cells.per.lineage.facet.condition <- function(){
  
  df <- read.table("01_metadata_exploration/03_Aquino2023-filtered-cell-count-01_lineage-Condition.txt")
  condition.vec <-  c("NS" = "Baseline", "IAV" = "IAV", "COV" = "COV")
  df$Condition <- plyr::revalue(df$Condition, condition.vec)
  df$Condition <- factor(df$Condition, levels = condition.vec)
  color.vec <- setNames(lineage.df$color, nm = lineage.df$lineage)
  df$lineage <- factor(df$lineage, levels = rev(names(color.vec)))
  # Params
  font_size = 6
  line_with = 0.25
  # Plot 
  x_lab <- "N cells"
  y_lab <- "Cell type"
  
  # df$label_color <- ifelse(df$lineage %in% c("T.CD4", "MONO"), "white", "black")
  # 
  n_cell.plt <- ggplot(data = df, aes(x = n_cells, y = lineage, color = lineage)) +
    geom_bar(stat= "identity", linewidth=0.25, alpha = 1, fill = "white") + 
    # Labels
    #geom_text(data = df, aes(x = n_cells, label = n_cells),
              # color = ifelse(lineage %in% c("T.CD4", "MONO"), yes = "white", no = "black")), # this doesnt work evet!
              #size = font_size / .pt, hjust = 0.5, position = position_stack(vjust = .5)) +
    facet_grid(Condition ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(breaks = c(25, 75) * 1e3, labels = c("25K", "75K")) + 
    scale_color_manual( values = lineage.color.vec ) +
    xlab(x_lab) + ylab(y_lab) + 
    plot_theme.aes + theme(axis.title.y = element_blank(), axis.title.x =element_text(size = 6))+
    guides(fill = FALSE, color = FALSE ) 
  
  n_cell.plt
}

n_cell.lineage.by.condition.plt <- generate.n.cells.per.lineage.facet.condition()

##  2.2 N of cells per Cell type  ------------------------------------------
generate.n.cells.per.celltype.facet.lineage <- function(){
  
  df <- read.table("01_metadata_exploration/03_Aquino2023-filtered-cell-count-02_celltype-Condition.txt")
  
  df$Condition <- plyr::revalue(df$Condition, condition.vec)
  df$Condition <- factor(df$Condition, levels = rev(condition.vec) )
  color.vec <- setNames(celltype.df$color, nm = celltype.df$celltype)
  df$celltype <- factor(df$celltype, levels = rev(unique(df$celltype)))
  df$lineage <- factor(df$lineage, levels = lineage.df$lineage)
  
  # Plot 
  x_lab <- "N cells"
  y_lab <- "Cell type"
  
  font_size = 6
  line_with = 0.25
  
  # 
  n_cell.plt <- ggplot(data = df, aes(x = n_cells, y = celltype, fill = Condition)) +
    geom_bar(stat= "identity", linewidth=0.25, alpha = 1) +  #, fill = "white") + 
    # Labels # disabling or now
    # geom_text(data = df, aes(x = n_cells, label = n_cells),
    #           # color = ifelse(cell_type %in% c("T.CD4", "MONO"), yes = "white", no = "black")), # this doesnt work evet!
    #           size = 3, hjust = 0.5, position = position_stack(vjust = .5)) +
    facet_grid(lineage ~ ., scales = "free_y", space = "free_y", switch = "y") +
     scale_x_continuous(breaks = c(50, 100) * 1e3, labels = c("50K", "100K")) + 
    scale_fill_manual( values = condition.color.vec ) +
     scale_color_manual( values = condition.color.vec ) +
    xlab(x_lab) + ylab(y_lab) + 
    theme(plot.title = element_blank(), plot.subtitle = element_blank(),
          panel.background = element_blank(), plot.background = element_blank(),
          axis.line.x = element_line(linewidth = line_with), axis.line.y = element_line(linewidth = line_with),
          axis.text.x = element_text(size = font_size, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = font_size), 
          axis.title.y = element_blank(),
          axis.title = element_text(size = font_size),
          strip.text = element_text(size = font_size),
          legend.position = c(0.9, 0.15), 
          # legend.position="bottom", legend.box = "vertical",
          legend.key.size = unit(0.1, "cm"), 
          legend.title=element_text(size=font_size-1), legend.text=element_text(size=font_size-2),
          strip.background = element_rect(fill = NA, color = "black", size = line_with), 
          strip.placement = "outside") + 
    guides( color = FALSE, fill = guide_legend(reverse=T) ) 
  
  n_cell.plt
}

n_cell.celltype.facet.lineage.plt <- generate.n.cells.per.celltype.facet.lineage()

generate.n.cells.per.celltype <- function(){
  
  df <- read.table("01_metadata_exploration/03_Aquino2023-filtered-cell-count-02_celltype-Condition.txt")
  # factor by overall celltype N of cells, not by cellytpe_Condition combination 
  count.df <- df %>% group_by(celltype) %>% summarise(n_cells = sum(n_cells)) %>% arrange(desc(n_cells)) 
  df$celltype <- factor(df$celltype, levels = rev(count.df$celltype))
  # subset to input cts 
  # if(!is.na(cts)){ df <- df[df$celltype %in% cts, ] }
  
  df$Condition <- plyr::revalue(df$Condition, condition.vec)
  df$Condition <- factor(df$Condition, levels = rev(condition.vec) )
  color.vec <- setNames(celltype.df$color, nm = celltype.df$celltype)
  df$lineage <- factor(df$lineage, levels = lineage.df$lineage)
  
  # Plot 
  x_lab <- "N cells"
  y_lab <- "Cell type"
  
  # 
  n_cell.plt <- ggplot(data = df, aes(x = n_cells, y = celltype, fill = Condition)) +
    geom_bar(stat= "identity", linewidth=0.25, alpha = 1, color = "white") +  #, fill = "white") + 
    # Labels # disabling or now
    # geom_text(data = df, aes(x = n_cells, label = n_cells),
    #           # color = ifelse(cell_type %in% c("T.CD4", "MONO"), yes = "white", no = "black")), # this doesnt work evet!
    #           size = 3, hjust = 0.5, position = position_stack(vjust = .5)) +
    # facet_grid(lineage ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(breaks = c(50, 100) * 1e3, labels = c("50K", "100K")) + 
    scale_fill_manual( values = condition.color.vec ) +
    scale_color_manual( values = condition.color.vec ) +
    xlab(x_lab) + ylab(y_lab) + 
    plot_theme.aes + theme(axis.title.y = element_blank(), axis.title.x = element_text(size =6),
                       legend.position = c(0.9, 0.15), legend.title = element_blank()) + 
    guides(fill = guide_legend(override.aes = list(label = "", size = 1, ncol=1), reverse = T))
  n_cell.plt
}

n_cell.celltype.plt <- generate.n.cells.per.celltype()
