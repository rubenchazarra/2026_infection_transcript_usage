#!/usr/bin/env Rscript 

# Over Representation Analysis (ORA) for differential transcript usage output ( Note: whether DTU genes are shared with DGE upregulated / downregulated is added previously in the 02_04_infectionDTU-processing.R script )
# ORA is calculated for different gene sets: i) all DTU genes, ii) DTU genes which are not DGE, iii) DTU genes which are DGE, iv) DTU genes which are DGE upregulated, v) DTU genes which are DGE downregulated

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--dtu"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'List of processed DTU output' 
  ), 
  make_option(
    c("-o", "--out_file_go"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file for ORA on GO terms' 
  ), 
  make_option(
    c("-u", "--out_file_msigdb"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file for ORA on Molecular Signature DataBase (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) ' 
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## FUNCTIONS ##
gene.set.creator.new <- function(df){
  list( "DTU.all" = df[ df$fdr < 0.05, "gene_id"],
        "DTU.unique" = df[ df$fdr < 0.05 & df$is.DGE == F, "gene_id"], 
        "DTU.DGE" = df[ df$fdr < 0.05 & df$is.DGE == T, "gene_id"], 
        "DTU.DGE.Up" = df[ df$fdr < 0.05 & df$is.DGE == T & df$DGE.direction == "Up", "gene_id"], 
        "DTU.DGE.Down" = df[ df$fdr < 0.05 & df$is.DGE == T & df$DGE.direction == "Down", "gene_id"]
  )
}


ora.clusterProfiler <- function( query, background, p.value.thres = 0.05 ){
  
  # remove ENSEMBL id version 
  query <- unique(gsub("\\..*", "", query))
  background <- unique(gsub("\\..*", "", background))
  
  # Iterate across databases
  db.vec <- c("BP")#, "MF", "CC") 
  
  l <- lapply(db.vec, function(db){
    print(paste0("Computing ORA for db: ", db))
    
    clusterProfiler::enrichGO(gene = query,
                              universe = background,
                              keyType = "ENSEMBL",
                              OrgDb = org.Hs.eg.db,
                              pvalueCutoff = p.value.thres,
                              pAdjustMethod = "BH", 
                              ont = db, 
                              minGSSize = 5,
                              readable = T)
  })
  names(l) <- db.vec
  return(l)
}

ora.clusterProfiler.wrapper <- function(df) {
  
  if(class(df) == "list") {lapply(df, ora.clusterProfiler.wrapper)
  }else{
    # 1) create gene sets
    all.tested.genes <- df$gene_id
    gene.set.l <- gene.set.creator.new(df)
    # 2) run ORA with DB = "BP"
    ora.l <- lapply(gene.set.l, function(gene.set) ora.clusterProfiler(query = gene.set, background = all.tested.genes, p.value.thres = 0.05))
    return(ora.l)
  }
}

msigDB.clusterProfiler <- function( query, background, msigdb.hs ) {
  print("Running msigDB clusterProfiler ...")
  # remove ensembl id version 
  query <- unique(gsub("\\..*", "", query))
  background <- unique(gsub("\\..*", "", background))
  
  clusterProfiler::enricher(query,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = background,
                            minGSSize = 5,
                            maxGSSize = 500,
                            qvalueCutoff = 0.2,
                            TERM2GENE = msigdb.hs,
                            TERM2NAME = NA
  )
}


msigDB.clusterProfiler.wrapper <- function(df, msigdb.hs) {
  
  if(class(df) == "list") { lapply(df, function(df) msigDB.clusterProfiler.wrapper(df, msigdb.hs))
  }else{
    # 1) create gene sets
    all.tested.genes <- df$gene_id
    gene.set.l <- gene.set.creator.new(df)
    # 2) run MSigDB enrichment
    ora.l <- lapply( gene.set.l, function(gene.set) msigDB.clusterProfiler(query = gene.set, background = all.tested.genes, msigdb.hs ))
    return(ora.l)
  }
}

## PACKAGES ##
suppressPackageStartupMessages(require(clusterProfiler))
suppressPackageStartupMessages(require(WebGestaltR))
suppressPackageStartupMessages(require(org.Hs.eg.db))
suppressPackageStartupMessages(require(msigdbr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(msigdb))
suppressPackageStartupMessages(require(msigdbr))
## PARAMS
dtu <- opt$dtu
out_file_go <- opt$out_file_go
out_file_msigdb <- opt$out_file_msigdb

###############
##### RUN ##### 
###############

dtu <- readRDS(dtu)
# select lineage 
dtu <- dtu$lineage

print("Input files read ...")

if(!is.na(out_file_go)){
  print("Running ORA GP:BP...")
  # ORA GO:BP
  ora.l <- ora.clusterProfiler.wrapper(dtu)
  print("ORA GP:BP computed ...")
  # Save
  saveRDS(ora.l, out_file_go)
  print("Outfile GP:BP saved ...")
}

if(!is.na(out_file_msigdb)){
  print("Running ORA MSigDB ...")
  # ORA MSigDB
  ## We can retrieve all human gene sets:
  msigdb.hs <- msigdbr::msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, ensembl_gene) # previously only doing in some categories
  # hallmark.db.l = list("H" = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, ensembl_gene), # H --> Hallmark gene sets, 50 gene sets
  #                      "C7" = msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, ensembl_gene) # C7 --> immunologic gene sets, 5219 gene sets )
  
  msig.db.l <- msigDB.clusterProfiler.wrapper(dtu, msigdb.hs)
  print("ORA MSigDB computed ...")
  # Save
  saveRDS(msig.db.l, out_file_msigdb)
  print("Outfile MSigDB saved ...")
}