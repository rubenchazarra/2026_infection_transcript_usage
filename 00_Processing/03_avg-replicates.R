#!/usr/bin/env RScript

# The Aquino et al.,(2023), contains technical replicates. This script averages counts from replicate samples of same Donor.ID and Condition
# In the origina publication they averaged the batch-corrected CPM values obtained across different replicates for the same individual and set of stimulation conditions, to obtain final estimates of gene expression

suppressPackageStartupMessages(require("optparse"))
option_list = list(
  make_option(
    c("-i", "--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Pseudobulk counts'
  ),
  make_option(
    c("-o", "--out_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output CPM normalized file suffixes'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
suppressPackageStartupMessages(require(dplyr))

### FUNCTIONS ###
avg_replicates <- function(counts, meta_data){
  # Average replicates from same Donor.ID based on sample_metadata, based on a Sample_ID composed of Library_POP_Donor.ID_Condition created in "../../00_MetaData/02_metadata_expansion.Rmd"
  # All Donor.ID_Condition replicates are in different libraries of the same Run, except for Run==16 which contains additional replicates used for "Quantification of batch effects and replicability" in Anquino 2023
  
  if(class(counts) == "list"){lapply(counts, function(counts) avg_replicates(counts, meta_data))
    
    }else{
    counts  = as.data.frame(t(counts))
    # IMPORTANT NOTE (2024-08-21) in principle Donor.ID-Condition replicates are in the same Run, except for Run16 in which we have Donor.IDs also present in another Runs, this is why we include Run, to avoid averaging Donors from different Runs
    counts = cbind( meta_data[ match(rownames(counts), meta_data$Sample.ID), c("Run", "Donor.ID", "Condition") ], counts)
    # Remove Run==16 data
    # IMPORTANT NOTE (2024-08-21) Excluding samples from Run16. We could also keep them and include Donor.ID in the DTU model by MANTA but MANTA does not handdle well factors with many levels such as Donor.ID
    # IMPORTANT NOTE (2024-09-13) After chatting with Aquino (author) I found Run==16 are replicates used to assess the replication of samples section ("Quantification of batch effects and replicability"), hence it is correct to exclude them 
    # print(" IMPORTANT NOTE (2024-08-21) Excluding samples from Run16. We could also keep them and include Donor.ID in the DTU model by MANTA but MANTA does not handdle well factors with many levels such as Donor.ID")
    # NOTE 2025-05-04 We are already considering the final samples: i) in_final_dataset, ii) Run != 16, iii) TP == "T6", iv) POP %in% c("EUB", "AFB"),
    # counts <- counts [ counts[["Run"]] != "16", ]
    # Average
    avg.counts <- as.data.frame( counts %>% group_by( Run, Donor.ID, Condition ) %>% summarise(across(everything(), mean)) ) # grouping by Run, Donor.ID, Condition because All Donor.ID_Condition replicates are in different libraries of the same Run
    print("Averaged counts across Donor.ID, Condition ...")
    
    # Rename samples as Donor.ID_Condition
    rownames( avg.counts ) <- paste( avg.counts[["Donor.ID"]], avg.counts[["Condition"]], sep = "_")
    
    # Remove metadata columns from counts
    avg.counts <- avg.counts [, -c ( match( c("Run", "Donor.ID", "Condition"),  colnames(avg.counts))) ]
    
    # Flip to columns = samples; rows = features
    avg.counts <- as.data.frame(t(avg.counts))
    
    return(avg.counts)
    }
}

# 0. Params
counts <- opt$counts
out_file <- opt$out_file

# 1. Read pseudobulk file
counts <- readRDS( counts )
meta_data <- read.table("../../00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)

# 2. Average replicates
avg.counts <- avg_replicates (counts, meta_data)

# 3. Save object
saveRDS(avg.counts, out_file)