#!/usr/bin/env Rscript

##  0. Purpose ------------------------------------------

## @RCG 2025-03-06 This script performs the initial processing of the infection DTU output data
# Please NOTE many hard coded output paths, so if changing input data, do change the output paths too to avoid overwriting of results

##  1. Functions ------------------------------------------
source("../../../../../../01_Randolph2021/05_Differential-Analyses/00_Scripts/00_DA-Functions.R")
source("../../../../../00_MetaData/03_Aesthetics-repo.R") # this load lineage and cell type df to link cell type to lineage to establish popDTU cell type specific
setwd("../03_Differential/02_DTU/02_DTU/01_infectionDTU/02_DA/")

##  2. Load data ------------------------------------------
file.l <- list("lineage" = read.RDS.files(path = "../01_infectionDTU/06_celltype-fractions/01_lineage/", pattern = "infectionDTU.rds"), 
               "celltype" = read.RDS.files(path = "../01_infectionDTU/05_Filt-1-5-0.8-0.01/02_celltype/", pattern = "infectionDTU.rds"))

##  2. Process ------------------------------------------

###  2.1 Select coef of interest ------------------------------------------
dtu.l <- select.df.from.nested.list(file.l, "Intercept")

###  2.2 Annotate ------------------------------------------

###  2.3 Add DGE info ------------------------------------------
dge.l <- readRDS("../../../../01_DGE/01_infectionDGE/02_DA/01_Processed-Data/05_include-Library/01_infectionDGE-filtered-list.rds")
# Ensure common lineages in both analysis
common.lineages <- intersect(names(dtu.l$lineage), names(dge.l$lineage))
dtu.l$lineage <- dtu.l$lineage[common.lineages]
dge.l$lineage <- dge.l$lineage[common.lineages]
# Ensure common lineages in both analysis
common.cell.types <- intersect(names(dtu.l$celltype), names(dge.l$celltype))
dtu.l$celltype <- dtu.l$celltype[common.cell.types]
dge.l$celltype <- dge.l$celltype[common.cell.types]

add.dge.info <- function(dtu, dge){
  if( class(dtu) == "list" & class(dge) == "list" ){ mapply(add.dge.info, dtu, dge, SIMPLIFY = F) # if list, iterate
  }else{
    # sharing
    dtu[["is.DGE"]] <- ifelse(dtu$gene_id %in% dge$genes, yes = T, no = F)
    # directionality
    dtu[["DGE.direction"]] <- ifelse( dge[ match(dtu$gene_id, dge$genes), "betas"]  > 0 , yes = "Up", no = "Down" )
    return(dtu)
  }
}
dtu.l <- add.dge.info(dtu.l, dge.l)

###  2.4 Add DTE info ------------------------------------------
# Finally not interested

###  2.5 Add Transcript switch info  ------------------------------------------
add.tx.switch.info <- function(dtu){
  print("Adding transcript switch information ...")
  if( class(dtu) == "list" ){ lapply(dtu, add.tx.switch.info) # if list, iterate
  }else{
    is.transcript.switch <- apply(dtu, 1, function(x){ # This is to compute if most abundant ECs share transcripts among them, or covnersly we can talk about real isoform swithcing
      x1 = x[["max.tx_name.1"]];  x2 = x[["max.tx_name.2"]]
      ec1.in.ec2 <- grep(x1, x2) ;  ec1.in.ec2 <- if(length (ec1.in.ec2)==0){F}else{T}
      ec2.in.ec1 <- grep(x2, x1); ec2.in.ec1 <- if(length (ec2.in.ec1)==0){F}else{T}
      # if EC1 is not in EC2 or viceversa, that is an equivalence class switch but not a transcrpt switch
      transcript.switch <- !(ec1.in.ec2 | ec2.in.ec1)
      return(transcript.switch)
    })
    dtu[["is.transcript.switch"]] <- unlist(is.transcript.switch)
    return(dtu)
  }
}
dtu.l <- add.tx.switch.info(dtu.l)  

###  2.5 Save as list  ------------------------------------------
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-unfiltered-list.rds"
saveRDS(dtu.l, file.name)

###  2.6 Save as DF  -----------------------------------------
# merge Condition
dtu.l <- lapply(dtu.l, function(l) lapply(l, function(ct)  plyr::ldply(ct, .id = "Condition")))
# merge celltype
dtu.l <- lapply(dtu.l, function(l) plyr::ldply(l, .id = "celltype"))
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-unfiltered-df.rds"
saveRDS(dtu.l, file.name)

###  2.7 Cell type specificity (compared to lineage)  BUG here -----------------------------------------
add.lineage.celltype.specificity <- function(file.name){
  print("1. Adding lineage specificity ...")
  dtu.l <- readRDS(file.name)
  lineage <- dtu.l$lineage
  lineage$lineage <- lineage$celltype # add lineage col to lineage
  celltype <- dtu.l$celltype
  lineage$is.DTU.lineage <- ifelse(lineage$fdr < 0.05, T, F)
  celltype$is.DTU.celltype <- ifelse(celltype$fdr < 0.05, T, F)
  celltype$lineage <- celltype.df[ match(celltype$celltype, celltype.df$celltype), "lineage"] # add lineage to celltype
  celltype <- celltype[celltype$is.DTU.celltype, ]
  ## 1) Add cell type specificity in lineage
  cols.celltype <- c("lineage", "celltype", "Condition", "gene_id", "gene_name", "is.DTU.celltype")
  celltype <- celltype[ , cols.celltype ] %>% group_by(lineage, Condition, gene_id, gene_name, is.DTU.celltype) %>% summarise(is.DTU.celltype.which = paste0(unique(celltype), collapse =","))
  lineage <- merge(lineage, celltype, by = c("lineage", "Condition", "gene_id", "gene_name"), all.x = T, all.y = F )
  lineage$is.DTU.lineage.celltype <- NA
  ## add specificity (only considering significant in lineage)
  lineage[ lineage$is.DTU.lineage  == T & is.na(lineage$is.DTU.celltype), "is.DTU.lineage.celltype"] <- "lineage.specific"
  lineage[ lineage$is.DTU.lineage  == T & !is.na(lineage$is.DTU.celltype) & lineage$is.DTU.celltype == T, "is.DTU.lineage.celltype"] <- "lineage.celltype"
  lineage.updated <- lineage
  
  ## 2) Add cell type specificity in celltype --> Here we want to know which genes are: celltype.specific, celltype.lineage shared and not.tested.in.lineage.
  print("Adding cell type specificity ...")
  lineage <- dtu.l$lineage
  celltype <- dtu.l$celltype
  lineage$is.DTU.lineage <- ifelse(lineage$fdr < 0.05, T, F)
  celltype$is.DTU.celltype <- ifelse(celltype$fdr < 0.05, T, F)
  celltype$lineage <- celltype.df[ match( celltype$celltype, celltype.df$celltype ), "lineage"] 
  # lineage <- lineage[lineage$is.DTU.lineage, ] # this uncomment is to identify "not.tested.in.lineage" category
  lineage$lineage <- lineage$celltype; lineage$celltype <- NULL # add lineage col to lineage and remove cell type (present in celltype)
  lineage$is.DTU.lineage.which <- lineage$lineage
  cols.lineage <- c("lineage", "Condition", "gene_id", "gene_name", "is.DTU.lineage")
  lineage <- lineage[ , cols.lineage ]
  celltype <- merge(celltype, lineage, by = c("lineage", "Condition", "gene_id", "gene_name"), all.x = T, all.y = F )
  celltype$is.DTU.lineage.celltype <- NA
  ## add specificity (only considering significant in lineage)
  celltype[ celltype$is.DTU.celltype == T & is.na(celltype$is.DTU.lineage), "is.DTU.lineage.celltype"] <- "not.tested.lineage"
  celltype[ celltype$is.DTU.celltype == T & ! is.na(celltype$is.DTU.lineage) & celltype$is.DTU.lineage == T, "is.DTU.lineage.celltype"] <- "lineage.celltype"
  celltype[ celltype$is.DTU.celltype == T & ! is.na(celltype$is.DTU.lineage) & celltype$is.DTU.lineage == F, "is.DTU.lineage.celltype"] <- "celltype.specific"
  celltype[ celltype$is.DTU.celltype == F & ! is.na(celltype$is.DTU.lineage) & celltype$is.DTU.lineage == F, "is.DTU.lineage.celltype"] <- "lineage.specific"
  celltype.updated <- celltype
  
  dtu.l$lineage <- lineage.updated; dtu.l$celltype <- celltype.updated
  print("Saving to input file name ...")
  saveRDS(dtu.l, file.name)
}

file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-unfiltered-df.rds"
add.lineage.celltype.specificity(file.name)

###  2.5 Save as list again (to include cell type specificity)  ------------------------------------------
## split by cell type
dtu.l <- readRDS(file.name)
dtu.l <- list("lineage" = split(dtu.l$lineage, dtu.l$lineage$lineage),
              "celltype" = split(dtu.l$celltype, dtu.l$celltype$celltype))
dtu.l <- lapply(dtu.l, function(x) lapply(x, function(ct) split(ct, ct$Condition )))
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-unfiltered-list.rds"
saveRDS(dtu.l, file.name)

###  2.6 Filter FDR  -----------------------------------------
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-unfiltered-list.rds"
file.l <- readRDS( file.name )
dtu.l <- filter.fdr(file.l, fdr.col = "fdr", fdr.thres = 0.05)

# Save as list --> for Visualization purposes
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-filtered-list.rds"
saveRDS(dtu.l, file.name)

###  2.7 Save as DF  -----------------------------------------
# merge Condition
dtu.l <- lapply(dtu.l, function(l) lapply(l, function(ct)  plyr::ldply(ct, .id = "Condition")))
# merge celltype
dtu.l <- lapply(dtu.l, function(l) plyr::ldply(l, .id = "celltype"))
file.name <- "01_Processed-Data/06_celltype-fractions/01_infectionDTU-filtered-df.rds"
saveRDS(dtu.l, file.name)

print("Processing of infectionDTU results finished !")
