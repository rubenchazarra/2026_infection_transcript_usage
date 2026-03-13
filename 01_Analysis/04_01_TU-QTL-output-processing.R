#!/usr/bin/env Rscript

##  0. Purpose ------------------------------------------

## This script performs initial processing of sQTLseekeR2 output for transcritp usage (TU) quantitative trait loci (QTL) output
# # Please NOTE many hard coded output paths, so if changing input data, do change the output paths too to avoid overwriting of results

# PACKAGES
require(dplyr)

# 0. 
tx2gene <- read.table("../../01_Quantif/00_IndexData/01_GENCODEv42-IAV-SARS_Cov_2-tx2gene-with-ENS.tsv", sep = "\t", header = T)
# 1. PERMUTED (SIGNIFICANT) results
cts <- c("B", "MONO", "NK", "T.CD4", "T.CD8")
conditions <- c("COV", "IAV", "NS")
path <- "../02_run-nf/02_results/01_lineage/02_celltype-fractions/"
ct.l <- lapply(cts, function(ct){
  file.path.ct <- paste0(path, ct)
  condition.l <- lapply(conditions, function(cond){
    file.path.ct.cond <- paste0(file.path.ct, "/03_sQTL_permuted/", cond, "/sqtls-0.01fdr.permuted.tsv")
    tab = read.table(file.path.ct.cond, sep = "\t", header = T)
    # add gene name
    tab$gene_name <- tx2gene[ match(tab$geneId, tx2gene$gene_id), "gene_name"]
    return(tab)
  })
  names(condition.l) <- conditions 
  return(condition.l)
})
names(ct.l) <- cts

## 1.1 Save as list 
saveRDS(ct.l, "01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-permuted-fdr-0.01-list.rds")
## 1.2 Save as dataframe
ct.l <- lapply(ct.l, function(ct) plyr::ldply(ct, .id="Condition"))
df <- plyr::ldply(ct.l, .id = "lineage")
saveRDS(df, "01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-permuted-fdr-0.01-df.rds")
print("Permuted results saved ...")

# 2. TESTED results
ct.l <- lapply(cts, function(ct){
  file.path.ct <- paste0(path, ct)
  condition.l <- lapply(conditions, function(cond){
    file.path.ct.cond <- paste0(file.path.ct, "/02_sQTL_nominal/", cond, "/all-tests.nominal.tsv")
    tab = read.table(file.path.ct.cond, sep = "\t", header = T)
    # add gene name
    tab$gene_name <- tx2gene[ match(tab$geneId, tx2gene$gene_id), "gene_name"]
    return(tab)
  })
  names(condition.l) <- conditions 
  return(condition.l)
})
names(ct.l) <- cts
## 2.1 Save as list 
saveRDS(ct.l, "01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-tested-list.rds")
## 2.2 Save as dataframe
ct.l <- lapply(ct.l, function(ct) plyr::ldply(ct, .id="Condition"))
df <- plyr::ldply(ct.l, .id = "lineage")
saveRDS(df, "01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-tested-df.rds")
print("Tested results saved ...")


# 3. Add info for ORA

## 3.1 
perm <- readRDS("01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-permuted-fdr-0.01-df.rds")
perm <- perm %>% distinct(lineage, Condition, geneId, gene_name)
perm$is.TU.QTL <- TRUE
## 3.2 Add DTU 
dtu <- readRDS("../../03_Differential/02_DTU/02_DTU/02_popDTU/02_DA/01_Processed-Data/05_Filt-1-5-0.8-0.01/01_popDTU-filtered-df.rds")
dtu <- dtu$lineage
dtu <- dtu %>% distinct(lineage, Condition, gene_name)
dtu$is.DTU <- TRUE
perm <- merge(perm, dtu, by = c("lineage", "Condition", "gene_name"), all.x = T, all.y = F)

tested <- readRDS("01_Processed-Data/03_celltype-fractions/01_sQTLseekeR2-tested-df.rds")
tested <- tested %>% distinct(lineage, Condition, geneId, gene_name)

tested <- merge(perm, tested, by = c("lineage", "Condition", "geneId", "gene_name"), all = T)
tested[is.na(tested)] <- F
tested.l <- split(tested, tested$lineage)
tested.l <- lapply(tested.l, function(x) split(x, x$Condition))
saveRDS(tested.l, "04_ORA/01_Processed-Data/03_celltype-fractions/00_sQTLseekeR-tested-genes.rds")

print("Finished collecting sQTLseekeR2 results...")