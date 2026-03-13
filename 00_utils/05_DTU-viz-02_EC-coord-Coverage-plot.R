#!/usr/bin/env Rscript

##  0. packages ------------------------------------------
suppressPackageStartupMessages(require(Gviz))

##  0. PALETTE ------------------------------------------

ec.palette.general <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#fca311", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", 
                         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#ffd670", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5")
other.ECs.col <- "#6c757d"

# 0.2 population DTU 
pop.color.vec <<- c( "EUB"="#212529", "AFB"="#9e0059" )

# 0.3 interaction DTU s
pop.condition.color.vec <<- c( "EUB_NS"="#6c757d", "EUB_IAV"="#212529", "EUB_COV"="#212529",
                              "AFB_NS"="#ff0054", "AFB_IAV"="#9e0059", "AFB_COV"="#9e0059" )

##  1. PROCESSING ------------------------------------------
pt_to_cex <- function(pt_size) { pt_size / 12 } # Convert point_size to cex value

# Set of functions for Equivalence class visualization
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

parse.gtf.gene.iso <- function(gtf, feature_name, set ){
  ## Parse GTF for visualization by Gviz
  if(set == "gene"){
    gtf.i <- gtf [ gtf[["gene_name"]] %in% feature_name, ]
  } else if(set == "transcript"){
    gtf.i <- gtf [ gtf[["transcript_name"]] %in% feature_name, ]
  }
  # Note: 'exon' type in GENCODE includes UTR+CDS. Non-coding features do not have UTR or CDS
  if(! all(c("CDS", "UTR") %in% gtf.i[["type"]] )) {
    transcript.model <- gtf.i[gtf.i$type %in% c("exon"), ]
  }else {
    transcript.model <- gtf.i[gtf.i$type %in% c("CDS","UTR"), ]
  }
  # Rename columns --
  transcript.model$gene <- transcript.model$gene_name
  transcript.model$exon <- transcript.model$exon_id
  transcript.model$transcript <- transcript.model$feature_name
  transcript.model$gene <- transcript.model$gene_id
  direction <- ifelse(transcript.model$strand == "+", "<", ">")
  transcript.model$symbol <- paste(transcript.model$feature_name, direction) # label in plot
  return(transcript.model)
}


##  2. TRANSCRIPT MODELS ------------------------------------------
tx.track.fun <- function(gen.coords, genome, feat.name){
  # Generate transcript track
  feat.type <- as.character(tolower(gen.coords$type)) # for different width of UTR and CDS, exon
  chr <- unique(as.character(gen.coords[["seqnames"]]))
  strand <- unique(as.character(gen.coords[["strand"]]))
  
  # Track
  track <- GeneRegionTrack(rstarts = gen.coords[["start"]],
                           rends = gen.coords[["end"]],
                           genome = genome,
                           chromosome = chr,
                           strand = strand,
                           name = feat.name, # this is the track name, appears on the left
                           transcriptAnnotation = "symbol",
                           symbol = feat.name, 
                           transcript = feat.name,
                           feature = feat.type, # distinguishes CDS from UTRs
                           fill = "#6c757d",     # 🔥 fill color for exons/arrows
                           col = "black",          # border color for exons/arrows
                           col.group = "black",    # color of transcript labels (names)
                           cex.group = pt_to_cex(10), 
                           showTitle = F, 
                           background.title = "transparent",  # makes the bar invisible
                           col.title = "transparent",          # also hide the text, 
                           fontfamily = "sans"  # same as ggplot2 default
                           )
  # Edit dpars (to enable distinction btw CDS and UTR)
  track@dp@pars[["collapse"]] <- FALSE
  return(track)
}

tx.track.generator <- function(gtf, gene_id, tx.annot.df ){
  
  tx_ids <- tx.annot.df[["transcript_id"]]
  tx_names <- tx.annot.df[["transcript_name"]]
  gtf.i <- gtf[gtf[["transcript_id"]] %in% tx_ids & gtf$type %in% c("exon", "UTR", "CDS"), ] # to later be able to visualize UTR as thinned box
  
  ## Generate list of transcript tracks
  # Transcript Tracks
  ## 1. GTF parsing
  gtf.l <- lapply( tx_names, function(tx_name) parse.gtf.gene.iso( gtf = gtf.i, feature_name = tx_name, set = "transcript" ))
  names(gtf.l) <- tx_names
  ## 2. Isoform Track generation
  tx.track.l <- lapply(tx_names, function(tx_name){ tx.track.fun(gen.coords = gtf.l[[tx_name]], genome = "hg38", feat.name = tx_name) } )
  names(tx.track.l) <- tx_names
  
  return(tx.track.l)
}

color.tx.track.fill <- function(tx.track.l, color = "black"){
  # add colors to tracks
  for( i in 1: length ( tx.track.l )) {
    tx.track.l[[i]]@dp@pars[["fill"]] <- color
    #tx.track.l[[i]]@dp@pars[["col"]] <- color
    tx.track.l[[i]]@dp@pars[["col.line"]] <- color
    #tx.track.l[[i]]@dp@pars[["fontcolor.group"]] <- color
  }
  return(tx.track.l)
}

##  3. EQUIVALENCE CLASS MODELS ------------------------------------------

EC.coordinate.generator <- function( gtf, gene_id, ecs ){
  #print("Extracting coordinates from equivalence classes ...")
  
  # All ECs
  ecs <- as.data.frame(ecs)
  ec.ids <- ecs [ ecs$gene_id %in% gene_id, "tcc" ]
  ec.names <- ecs [ ecs$gene_id %in% gene_id, "transcript_name" ]
  # list txs involved in each ECs
  ec.names.l <- setNames(lapply(ec.names, function(x) unlist(strsplit(x, ","))), nm = ec.names)
  
  # Create GTF
  gtf <- gtf[gtf$type %in% "exon", ] # keep only exons
  gtf.gr <- GRanges(seqnames = gtf$seqnames,
                    ranges = IRanges(start = gtf$start, end = gtf$end),
                    strand = gtf$strand,
                    transcript_id = gtf$transcript_id, 
                    transcript_name = gtf$transcript_name,
                    exon_number = gtf$exon_number
                    )

  # 1.  Common region of all Tx present in each TCC
  tx.intersect.l <- lapply(ec.names, function(ec){
    tx.vec <- unlist(strsplit(ec, ","))
    gtf.gr.i <- gtf.gr[gtf.gr$transcript_name %in% tx.vec , ]
    # Split GRanges by transcript ID
    gr.l <- split(gtf.gr.i, gtf.gr.i$transcript_name)
    
    # Find intersection transcripts
    common_regions <- Reduce(GenomicRanges::intersect, gr.l)
    return(common_regions)
  })
  names(tx.intersect.l) <- ec.names
  # 2. if and ECs is inside another one (such as tx_1,tx_2 is inside tx_1,tx_2,tx_3) substract
  ec.l <- lapply(ec.names, function(ec) {
    # in which other ECs are the tx of ec present ?
    ec.tx.l <- ec.names.l [[ ec ]]
    ec.tx.rest.l <- ec.names.l [ - match(ec, ec.names) ]
    presence.vec <- lapply(ec.tx.rest.l, function(x){ all(ec.tx.l %in% x) })
    presence.vec <- unlist(presence.vec)
    ec.query <- tx.intersect.l[[ec]]
    # if EC is in any other EC, proceed. In reality we are asking if the Tx in EC are among the Tx of any other EC
    if( any( presence.vec ) ){ 
      ec.rest <- tx.intersect.l[ match(names(presence.vec), names(tx.intersect.l)) ][ unname ( presence.vec ) ]
      # to generate the ec, iterativelly substract it to every other ECs with shared transcripts 
      for (i in ec.rest ) { ec.query <- GenomicRanges::setdiff(ec.query, i) }
    }
    # if output EC coordinates are empty, retrieve empty
    # if(length(ec.query) == 0){ ec.query <- GRanges( seqnames = NULL,
    #                                                 ranges = IRanges(start = NULL, end = NULL),
    #                                                 strand = NULL )
    # }
    return(ec.query)
    })
  names(ec.l) <- ec.names
  
  return(ec.l)  
}

EC.highlight.fun <- function(tracklist, ec.gtf, genome = "hg38"){
  
  ec.gtf <- as.data.frame(ec.gtf)
  chr <- unique(as.character(ec.gtf[["seqnames"]]))
  feat.name <- ec.gtf[["feature_name"]]
  
  ht <- HighlightTrack(trackList = unlist(tracklist, recursive = F), 
                       start = ec.gtf[["start"]], end = ec.gtf[["end"]],
                       chromosome = chr, 
                       name = feat.name, 
                       symbol = feat.name,
                       legend = TRUE,
                       cex.legend = pt_to_cex(1), # font size of the legend
                       fontfamily.legend = "sans", 
                       ncol = 1
                       )
  return(ht)
}

edit.EC.highlight.track.fill <- function(track.l, ec.annot.df){
  #print("Editting Equivalence class highlight track colors ...")
  
  track.names <- track.l@name
  track.name.df <- data.frame(track.name = track.l@name)

  # assign colors to tested ECs
  ec.palette <- ec.palette.general
  other.ECs.col <- other.ECs.col
  
  # add colors while dealing with lowly abundant ECs which are labelled as "other.ECs" in 01_boxplots and have a different color
  color.df <- ec.annot.df[, c("transcript_name", "exclude_plot")]
  color.df$color <- NA
  include.ecs <- color.df[ !color.df$exclude_plot, "transcript_name" ]
  color.df[ ! color.df$exclude_plot, "color" ] <- ec.palette[1: length(include.ecs)] 
  color.df[ color.df$exclude_plot, "color" ] <- other.ECs.col
  # merge both
  color.df <- merge(track.name.df, color.df, by.x = "track.name", by.y = "transcript_name", all = T) 
  # merge alters track order # 2025-02-17 tested.ecs for which we found no coordinates will be removed at this step. This is an issue we have to address
  color.df <- color.df[ match(track.name.df$track.name, color.df$track.name), ]
  # grey color for the non tested
  color.df[ is.na(color.df$exclude_plot), "exclude_plot"] <- "not.tested.ECs"
  color.df[ is.na(color.df$color), "color"] <- "#adb5bd"
  
  # add colors to tracks
  track.l@dp@pars[["fill"]] <-  color.df$color
  track.l@dp@pars[["col"]] <-  color.df$color
  # reduce alpha for the other tracks to be visible
  track.l@dp@pars[["alpha"]] <- 0.3 #0.5
  track.l@dp@pars$alpha.title = 0.3 #0.5
  # place EC highlight box in backgroud
  track.l@dp@pars$inBackground <- TRUE
  return(track.l)
}

##  4. Coverage Track Functions ------------------------------------------
avg.cov.of.replicates <- function(cov, group.vec ){
  # Example grouping vector (same length as number of columns)
  # Unique groups
  groups <- unique( group.vec )
  # Initialize a new dataframe to store the means
  result <- data.frame(matrix(ncol = length(groups), nrow = nrow(cov)))
  colnames(result) <- groups
  
  # Calculate the row-wise mean for each group
  for (group in groups) {
    result[[group]] <- rowMeans(cov[, group.vec == group, drop = FALSE])
  }
  # View result
  return(result)
}

aggregate.by.group.vec <- function(cov, group.vec, type){
  # Aggregate coverage file by SUM or MEAN
  if(!type %in% c("sum", "mean")) stop("type argument of aggregate.by.group.vec function must be one of: 'sum' or 'mean' ")
  # Unique groups
  groups <- unique ( group.vec )
  # Initialize a new dataframe to store the means
  result <- data.frame(matrix(ncol = length(groups), nrow = nrow(cov)))
  colnames(result) <- groups
  
  # Calculate the row-wise mean for each group
  for (group in groups) {
    if(type == "sum"){ result[[group]] <- rowSums(cov[, group.vec == group, drop = FALSE]) }
    if(type == "mean"){  result[[group]] <- rowMeans(cov[, group.vec == group, drop = FALSE]) }
  }
  # View result
  return(result)
}

read.cov.file <- function(path, celltype, gene_name, gene_id){
  # Read coverage file (previously computed with samtools depth)
  cov.file.name <- paste0(path, celltype, "_", gene_name, "-", gene_id, "-samtools-depth.txt.gz")
  if(! file.exists(cov.file.name)){ stop (paste0("File does not exist ....: ", ct, " - ", gene_name, " - ", gene_id, " ..."))
  }else{
    cov <- read.table(cov.file.name, sep = "\t", header = F, stringsAsFactors = FALSE)
    # add colnames to file
    bam.list.file.path <- paste0("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/02_cellranger/01_GENCODEv42-IAV-COV/03_coverage-bams/03_split-bams/02_bam-list-files/01_lineage/",
                                 celltype, "-cram-file-list.txt")
    bam.list.file.path <- read.table(bam.list.file.path, sep = "\t", header = F)
    # edit cov file names
    colnames(cov) <- c("chr", "pos", gsub("*.cram", "", basename(bam.list.file.path$V1)))
    colnames(cov) <- gsub(paste0("_", celltype), "", names(cov)) # remove cell type from naming to later operate with library_Donor.ID_Condition IDs
    return(cov) 
  }
}

bin_coverage_median <- function(df, bin_size = 10) {
  library(dplyr)
  # Identify columns to sum (exclude chr, start, end)
  cov_cols <- setdiff(colnames(df), c("chr", "pos"))
  # Compute bin start and sum all coverage columns
  df %>%
    mutate(pos = floor(pos/ bin_size) * bin_size) %>%
    group_by(chr, pos) %>%
    summarise(across(all_of(cov_cols), median), .groups = "drop")
}

process.cov.infectionDTU <- function(cov, condition.levels, condition.col, min.limit, max.limit, type){
  #print("Consider merging this function with the previous one ...")
  # 1. Select samples by condition
  annot.cols <- c("chr", "pos")
  annot.cov <- cov[, annot.cols]
  cov <- cov[, grep(paste0(condition.levels, collapse = "|"), colnames(cov)) ] # retain depth values for condition.levels only
  # 2. Average replicates
  sample_meta <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
  #sample_meta <- read.table("../../../../../00_MetaData/sample_metadata_expanded.tsv")
  sample_meta[["Library_Donor.ID_Condition"]] <- paste0(sample_meta$Library, "_", sample_meta$Donor.ID, "_", sample_meta$Condition)
  sample_meta.i <- sample_meta [ match(colnames(cov), sample_meta[["Library_Donor.ID_Condition"]]), ]
  cov <- avg.cov.of.replicates(cov, sample_meta.i [["Donor.ID_Condition"]] )
  ## Update metadata
  sample_meta.i <- sample_meta.i [ match(colnames(cov), sample_meta.i[["Donor.ID_Condition"]]), ]
  ## 3. Sum coverage by groups
  cov <- aggregate.by.group.vec(cov, sample_meta.i[[condition.col]], type = type)
  cov <- cbind( data.frame("chr" = annot.cov$chr, "start" = annot.cov$pos, "end" = annot.cov$pos + 1), cov)
  ## 4. Restric Cov to the coordinates
  cov <- cov[cov$start > min.limit & cov$end < max.limit, ]
  # 5. TODO: Does removing the zeroes affect the plot ? Objects will be lighter
  return(cov)
}

process.cov.interactionDTU <- function(cov, condition.levels, condition.col, min.limit, max.limit, type){
  # This function: i) averages replicated samples, ii) computer mean sample coverage for EUB_Condition_1, EUB_Condition_2, AFB_Condition_1, AFB_Condition_2
  annot.cols <- c("chr", "pos")
  annot.cov <- cov[, annot.cols]
  cov <- cov[, grep(paste0(condition.levels, collapse = "|"), colnames(cov)) ] # retain depth values for condition.levels only
  # 2. Average replicates
  sample_meta <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/00_MetaData/sample_metadata_expanded.tsv", sep = "\t", header = T)
  #sample_meta <- read.table("../../../../../00_MetaData/sample_metadata_expanded.tsv")
  sample_meta[["Library_Donor.ID_Condition"]] <- paste0(sample_meta$Library, "_", sample_meta$Donor.ID, "_", sample_meta$Condition)
  sample_meta.i <- sample_meta [ match(colnames(cov), sample_meta[["Library_Donor.ID_Condition"]]), ]
  cov <- avg.cov.of.replicates(cov, sample_meta.i [["Donor.ID_Condition"]] )
  ## Update metadata
  sample_meta.i <- sample_meta.i [ match(colnames(cov), sample_meta.i[["Donor.ID_Condition"]]), ]
  sample_meta.i[["POP_Condition"]] <- factor(paste0(sample_meta.i$POP, "_", sample_meta.i$Condition ), levels = c("AFB_NS", "AFB_IAV", "AFB_COV", "EUB_NS", "EUB_IAV", "EUB_COV"))
  ## 3. Sum coverage by groups
  # condition.col <- "POP_Condition"
  cov <- aggregate.by.group.vec(cov, sample_meta.i[[condition.col]], type = type)
  cov <- cbind( data.frame("chr" = annot.cov$chr, "start" = annot.cov$pos, "end" = annot.cov$pos + 1), cov)
  ## 4. Restric Cov to the coordinates
  cov <- cov[cov$start > min.limit & cov$end < max.limit, ]
  return(cov)
}

coverage.track.fun <- function(cov, groups, groups.color, track.name){
  ## Generate a single Gviz data track with coverage for groups and coloured by groups.color
  # Create GRanges object
  gr <- GRanges(cov)
  # 2025-07-29 This does not work ! I'm reordeing the cov dataframe from before
  # NOTE: groups are always ordered in alphabetical order in plot legend independently of order specified in groups argument
  # We have to harmonize the groups and the groups color so that they don't get flipped
  # Set groups as factor with levels in desired order
  # groups <- factor(groups, levels = groups)
  # Reorder the color vector accordingly
  # groups.color <- groups.color[ groups ]
  # DataTrack
  dataTrack <- DataTrack(range = gr, 
                         genome = "hg38", 
                         type = "l",  
                         name = track.name, 
                         groups = groups, 
                         col = groups.color, 
                         lwd = 1, 
                         col.axis="black", col.title="black", # axis color
                         cex.axis = pt_to_cex(12), cex.title = pt_to_cex(14), # font size
                         fontface.title = 1,   # 1 = plain, 2 = bold, 3 = italic # font aspect
                         fontface.axis = 1, 
                         background.title="transparent",  # axis bar
                         legend = T,
                         fontcolor.legend = "black",   # legend
                         fontface.legend = 1, 
                         fontfamily.legend = "sans", 
                         fontsize.legend=14, 
                         box.legend = TRUE
                         )
  return(dataTrack)
}

##  5. Wrapper functions  ------------------------------------------
infectionDTU.Tx.Cov.EqClass.track.wrapper <- function( dtu.l, rel.l, ct.condition ){
  
  # 0. Get DTU and Count data
  dtu.df = dtu.l[[ ct.condition ]] 
  dtu.df <- dtu.df[dtu.df$n_features < 20, ] # fitler, we don't want to plot these
  rel_exp = rel.l[[ ct.condition ]] # Dataframe including 2 conditions, no longer a list of conditions
  # Remove annotation from diff_exp and rel_exp (we will reannotate later)
  rel_exp <- rel_exp[, intersect(colnames(rel_exp), meta_data[["Donor.ID_Condition"]])]
  
  ct <- gsub(".NS|.COV|.IAV", "", ct.condition)
  condition.pairwise <- gsub(paste0(ct,"."), "", ct.condition) # if we modify the labels we fuck up downstream functions
  condition_1 <- strsplit(condition.pairwise, "[.]")[[1]][1] 
  condition_2 <- strsplit(condition.pairwise, "[.]")[[1]][2]
  
  track.l <- lapply( dtu.df[["gene_id"]], function(gene_id) { # iterate across genes
    
    tryCatch({ # To deal with errors
      
    i <- match( gene_id, dtu.df[["gene_id"]] )
    # 1. Annotate EC ids
    ec.vec <- unlist(strsplit(dtu.df[i , "feature_ids" ], ";"))
    if( ! dtu.df[i, "feature_ids.excluded"] == "NA" ) {  ec.vec <- c( ec.vec, unlist(strsplit(dtu.df[i , "feature_ids.excluded" ], ";"))) }
    ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 112 )
    gene_name <- unique(ec.annot.df$gene_name) 
    ec.annot.df <- order.ec.df( ec.annot.df )
    gene_name <- unique(ec.annot.df$gene_name) 
    ecs <- ec.annot.df$ec
    print(paste0(ct.condition, " - ", i, " - ", gene_name, " - ", gene_id))
    
    # 2. Subset counts
    ## 2.1 Annotate and Subset Rel Counts
    rel_exp.i <- rel_exp[ecs, ] # subset
    rel_exp.i <- cbind(annotate.equivalence.classess(ecs, tx2gene, n_cores = 4), rel_exp.i)
    
    ## 1.1 Collect lowly abundant ECs to "other.ECs" label (to match colouring in 01_Boxplot)
    annot.cols <- c("gene_id", "gene_name", "transcript_name", "ec")
    ec.mean.abundances <- data.frame("EC.mean.abundance" = rowMeans(rel_exp.i[, intersect(meta_data$Donor.ID_Condition, colnames(rel_exp.i)) ], na.rm = T))
    ec.annot.df <- cbind(ec.annot.df, ec.mean.abundances [ ec.annot.df$ec, , drop = F ] )
    if( nrow(ec.annot.df) > 6 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.abundance < 0.05
    } else { ec.annot.df$exclude_plot <- F }
    ## 3.2 If lowly abundant EC present, collapse into "other.ECs" label # 2025-07-29 I think we don't need this, we have already labeled the other.ECs in the ec.annot.df
    # if( length( ec.annot.df$exclude_plot[ec.annot.df$exclude_plot] ) > 1 ){ # only for rel_abundances 
    #   raw_exp.i <- out.l[["raw_exp.i"]] ; rel_exp.i <- out.l[["rel_exp.i"]] ; diff_exp.i <- out.l[["diff_exp.i"]] ; ec.annot.df <- out.l[["ec.annot.df"]]
    # }
    
    
    # 2. Tx Tracks
    ## 2.0 Annotation
    tx.vec <- unique(unlist(strsplit(ec.annot.df$ec, ",")))
    tx.annot.df <- annotate.equivalence.classess( tx.vec, tx2gene )
    names(tx.annot.df)[ grep("ec", names(tx.annot.df)) ] <- "transcript_id"
    tx.annot.df <- tx.annot.df[ order(tx.annot.df$transcript_name), ]
    ## 2.1. Generate Tx tracks
    tx.track.l <- tx.track.generator(gtf, gene_id, tx.annot.df) 
    ## 2.2. Edit tx.track colors --> Now editing inside GeneRegioTrack function
    # tx.track.l <- color.tx.track.fill(tx.track.l, color = "#6c757d")
    
    # 3. Coverage Tracks
    ## 3.1 Read coverage
    path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/02_cellranger/01_GENCODEv42-IAV-COV/03_coverage-bams/04_samtools-depth/02_samtools-depth/01_lineage/"
    cov <- read.cov.file( path, ct, gene_name, gene_id )
    # Median coverage across bins, this way the plot is better
    cov <- bin_coverage_median(cov, bin_size = 10)
    ## 3.2 Process coverage
    min.limit <- min( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "start"] )
    max.limit <- max( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "end"] )
    cov <- process.cov.infectionDTU(cov, 
                                    condition.levels = c(condition_1, condition_2), 
                                    condition.col = "Condition", 
                                    min.limit,
                                    max.limit, 
                                    type = "mean")
    ## 3.3 Coverage Track
    condition.color.vec <- c( "NS" = "#969696", "IAV" = "#75A9D1", "COV" = "#c1121f" )
    ### Reorder cols
    groups <- c(condition_1, condition_2)
    cov <- cov[, c("chr", "start", "end", groups)]
    cov.track <-  coverage.track.fun(cov,
                                     groups = c(condition_1, condition_2), 
                                     groups.color = c(condition.color.vec[condition_1], condition.color.vec[condition_2]), 
                                     track.name = "Mean coverage")
    
    # 4. Merge
    track.l <- c(tx.track.l, cov.track )
    
    # 5. Equivalence class Track
    ## 5.1 get EC coordinates (from the entire EC set)
    ec.coord.l <- EC.coordinate.generator( gtf, gene_id, ecs = all.ecs )
    ec.names <- names(ec.coord.l)
    ec.coord.l <- lapply( ec.names, function (ec) {
      ec.coord <- ec.coord.l[[ec]]
      if(length(ec.coord) > 0){  ec.coord$feature_name <- ec } # add feature name to GRanges
      return(ec.coord)
    })
    names(ec.coord.l) <- ec.names
    ec.coord.gtf <- Reduce(c, ec.coord.l)
    
    ## 5.2 ECs tested in DTU will be coloured and those not tested will be in grey or some other color outside of the pallette
    tested.ecs <- ec.annot.df$transcript_name
    
    ## 5.3 highlight tracks must be overlayed on top of an existing track list (tx.track and coverage track in our case)
    track.l <- EC.highlight.fun( tracklist = track.l, ec.coord.gtf )
    ## 5.4 Color tracks
    track.l <- edit.EC.highlight.track.fill( track.l, ec.annot.df )
    # close tryCatch
    }, error = function(e) e) 
  })  
  # add names
  gene_name.vec <- tx2gene[ match(dtu.df[["gene_id"]], tx2gene$gene_id), "gene_name"]
  names(track.l) <- gene_name.vec
  return(track.l)
}

popDTU.Tx.Cov.EqClass.track.wrapper <- function( dtu.l, ct.condition ){
  
  # 0. Get DTU and Count data
  dtu.df <- dtu.l[[ ct.condition ]] # ; dtu.df[is.na(dtu.df)] <- F # 2025-02-19 BUG, DTU genes not appearing in the DGE or DTE shoudl appear  as F, not as NA, edit in DA script
  dtu.df <- dtu.df[dtu.df$n_features < 20, ] # fitler, we don't want to plot these
  # 2025-04-05 Do we care about DGE or DTE to generate this plot ? 
  # rel_exp.l = rel.l[[ ct.condition ]] # list: IAV, NS # 2025-04-05 Not considering rel exp for now
  
  # Get celltype and condition
  # ct <- strsplit(ct.condition, "[.]")[[1]][1] # gsub("_.*","", ct.condition)
  ct <- gsub(".NS|.IAV|.COV", "", ct.condition)
  # condition <- strsplit(ct.condition, "[.]")[[1]][2] # gsub("*._","", ct.condition)
  condition <- gsub(paste0( ct, "."), "", ct.condition)
  
  track.l <- lapply( dtu.df[["gene_id"]], function(gene_id) { # iterate across genes
    
    tryCatch({ # To deal with errors
    
    i <- match( gene_id, dtu.df[["gene_id"]] )
    # 1. Annotate EC ids
    ec.vec <- unlist(strsplit(dtu.df[i , "feature_ids" ], ";"))
    ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 12 )
    ec.annot.df <- order.ec.df( ec.annot.df )
    gene_name <- unique(ec.annot.df$gene_name) 
    print(paste0(ct.condition, " - ", i, " - ", gene_name, " - ", gene_id))
    
    ## 1.1 Collect lowly abundant ECs to "other.ECs" label (to match colouring in 01_Boxplot) # 2025-04-05 We dont do in the boxplot either.
    # annot.cols <- c("gene_id", "gene_name", "transcript_name", "tcc")
    # rel_exp.l <- lapply( rel_exp.l, function(x) x[ ec.annot.df$ec, ]) # subset
    # ec.mean.abundances <- lapply(rel_exp.l, function(x) data.frame("EC.mean.abundance" = rowMeans(x[, !colnames(x) %in% annot.cols], na.rm = T)))
    # ec.mean.abundances <- setNames( cbind( ec.mean.abundances[[1]], ec.mean.abundances[[2]] ), nm = c("EC.mean.1", "EC.mean.2") ) # should match because we are working with the intersection of ECs in condition_1 and condition_2
    # ec.annot.df <- cbind(ec.annot.df, ec.mean.abundances [ ec.annot.df$ec, ] )
    # if( nrow(ec.annot.df) > 6 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.1 < 0.025 & ec.annot.df$EC.mean.2 < 0.025 
    # } else { ec.annot.df$exclude_plot <- F }
    ec.annot.df$exclude_plot <- FALSE # 2025-04-05 This line must be present for edit.EC.highlight.track.fill()
    
    # 2. Tx Tracks
    ## 2.0 Annotation
    tx.vec <- unique(unlist(strsplit(ec.annot.df$ec, ",")))
    tx.annot.df <- annotate.equivalence.classess( tx.vec, tx2gene )
    names(tx.annot.df)[ grep("ec", names(tx.annot.df)) ] <- "transcript_id"
    tx.annot.df <- tx.annot.df[ order(tx.annot.df$transcript_name), ]
    ## 2.1. Generate Tx tracks
    tx.track.l <- tx.track.generator(gtf, gene_id, tx.annot.df) 
    ## 2.2. Edit tx.track colors --> Now editing inside GeneRegioTrack function
    # tx.track.l <- color.tx.track.fill(tx.track.l, color = "#6c757d")
    
    # 3. Coverage Tracks
    ## 3.1 Read coverage
    path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/02_cellranger/01_GENCODEv42-IAV-COV/03_coverage-bams/04_samtools-depth/02_samtools-depth/01_lineage/"
    cov <- read.cov.file( path, ct, gene_name, gene_id )
    # Median coverage across bins, this way the plot is better
    cov <- bin_coverage_median(cov, bin_size = 10)
    ## 3.2 Process coverage
    min.limit <- min( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "start"] )
    max.limit <- max( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "end"] )
    cov <- process.cov.infectionDTU(cov, 
                                    condition.levels = condition, 
                                    condition.col = "POP", 
                                    min.limit,
                                    max.limit, 
                                    type = "mean")
    
    ## 3.3 Coverage Track
    groups.i <- setdiff(colnames(cov), c("chr", "start", "end"))
    ordered.groups <- c("AFB", "EUB") 
    groups.i <- ordered.groups[ordered.groups %in% groups.i]
    cov <- cov[, c("chr", "start", "end", ordered.groups)] ### Reorder cols
    cov.track <-  coverage.track.fun(cov,
                                     groups = groups.i, 
                                     groups.color = pop.color.vec[ groups.i ], 
                                     track.name = "Mean coverage")
    
    
    # 4. Merge
    track.l <- c(tx.track.l, list("Mean cov" = cov.track ))
    
    # 5. Equivalence class Track
    ## 5.1 get EC coordinates (from the entire EC set)
    ec.coord.l <- EC.coordinate.generator( gtf, gene_id, ecs = all.ecs )
    ec.names <- names(ec.coord.l)
    ec.coord.l <- lapply( ec.names, function (ec) {
      ec.coord <- ec.coord.l[[ec]]
      if(length(ec.coord) > 0){ # add feature name to GRanges
        ec.coord$feature_name <- ec
      }
      return(ec.coord)
    })
    names(ec.coord.l) <- ec.names
    ec.coord.gtf <- Reduce(c, ec.coord.l)
    
    ## 5.2 ECs tested in DTU vs not tested --> No need to split, we differentiate by the color
    # tested.ecs <- intersect( names(ec.coord.l), ec.annot.df$transcript_name )
    # tested.ecs <- tested.ecs[ match( ec.annot.df$transcript_name, tested.ecs )]
    tested.ecs <- ec.annot.df$transcript_name
    
    ## 5.3 highlight tracks must be overlayed on top of an existing track list (tx.track and coverage track in our case)
    track.l <- EC.highlight.fun( tracklist = track.l, ec.coord.gtf )
    ## 5.4 Color tracks
    track.l <- edit.EC.highlight.track.fill( track.l, ec.annot.df )
    return(track.l) 
    # close tryCatch
    }, error = function(e) e) 
  })
  # add names
  gene_name.vec <- tx2gene[ match(dtu.df[["gene_id"]], tx2gene$gene_id), "gene_name"]
  names(track.l) <- gene_name.vec
  return(track.l)
}


interactionDTU.Tx.Cov.EqClass.track.wrapper <- function( dtu.l, ct.condition ){
  
  # 0. Get DTU and Count data
  dtu.df = dtu.l[[ ct.condition ]] ; dtu.df[is.na(dtu.df)] <- F # 2025-02-19 BUG, DTU genes not appearing in the DGE or DTE shoudl appear as F, not as NA, edit in DA script
  # Avoid plotting genes with too many ECs
  dtu.df <- dtu.df[dtu.df$n_features < 20, ]
  ct <- gsub(".IAV.NS|.COV.NS|.IAV.COV", "", ct.condition)
  condition.pairwise <- gsub(paste0(ct,"."), "", ct.condition)
  # condition.pairwise.vec.short <- c("IAV.NS" = "IAV-Baseline",  "COV.NS" = "COV-Baseline", "IAV.COV" = "IAV-COV") 
  # condition.pairwise <- plyr::revalue(condition.pairwise, condition.pairwise.vec.short)
  condition_1 <<- strsplit(condition.pairwise, "[.]")[[1]][1] # Set as global variables to be used in add.meta_data()
  condition_2 <<- strsplit(condition.pairwise, "[.]")[[1]][2]
  
  track.l <- lapply( dtu.df[["gene_id"]], function(gene_id) { # iterate across genes
    
    tryCatch({ # To deal with errors
      
    i <- match( gene_id, dtu.df[["gene_id"]] )
    # 1. Annotate EC ids
    ec.vec <- unlist(strsplit(dtu.df[i , "feature_ids" ], ";"))
    if( ! dtu.df[i, "feature_ids.excluded"] == "NA" ) {  ec.vec <- c( ec.vec, unlist(strsplit(dtu.df[i , "feature_ids.excluded" ], ";"))) }
    ec.annot.df <- annotate.equivalence.classess( ec.vec, tx2gene, n_cores = 112 )
    ec.annot.df <- order.ec.df( ec.annot.df )
    gene_name <- unique(ec.annot.df$gene_name) 
    print(paste0(ct.condition, " - ", i, " - ", gene_name, " - ", gene_id))
    
    ## 1.1 Collect lowly abundant ECs to "other.ECs" label (to match colouring in 01_Boxplot) # 2025-04-30 Not doing this
    # annot.cols <- c("gene_id", "gene_name", "transcript_name", "tcc")
    # rel_exp.l <- lapply( rel_exp.l, function(x) x[ ec.annot.df$ec, ]) # subset
    # ec.mean.abundances <- lapply(rel_exp.l, function(x) data.frame("EC.mean.abundance" = rowMeans(x[, !colnames(x) %in% annot.cols], na.rm = T)))
    # ec.mean.abundances <- setNames( cbind( ec.mean.abundances[[1]], ec.mean.abundances[[2]] ), nm = c("EC.mean.1", "EC.mean.2") ) # should match because we are working with the intersection of ECs in condition_1 and condition_2
    # ec.annot.df <- cbind(ec.annot.df, ec.mean.abundances [ ec.annot.df$ec, ] )
    # if( nrow(ec.annot.df) > 6 ){ ec.annot.df$exclude_plot <- ec.annot.df$EC.mean.1 < 0.025 & ec.annot.df$EC.mean.2 < 0.025 
    # } else { ec.annot.df$exclude_plot <- F }
    ec.annot.df$exclude_plot <- F
    
    # 2. Tx Tracks
    ## 2.0 Annotation
    tx.vec <- unique(unlist(strsplit(ec.annot.df$ec, ",")))
    tx.annot.df <- annotate.equivalence.classess( tx.vec, tx2gene )
    names(tx.annot.df)[ grep("ec", names(tx.annot.df)) ] <- "transcript_id"
    tx.annot.df <- tx.annot.df[ order(tx.annot.df$transcript_name), ]
    ## 2.1. Generate Tx tracks
    tx.track.l <- tx.track.generator(gtf, gene_id, tx.annot.df) 
    ## 2.2. Edit tx.track colors --> Now editing inside GeneRegioTrack function
    # tx.track.l <- color.tx.track.fill(tx.track.l, color = "#6c757d")
    
    # 3. Coverage Tracks
    ## 3.1 Read coverage
    path = "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/01_Quantif/02_cellranger/01_GENCODEv42-IAV-COV/03_coverage-bams/04_samtools-depth/02_samtools-depth/01_lineage/"
    cov <- read.cov.file( path, ct, gene_name, gene_id )
    # Median coverage across bins, this way the plot is better
    cov <- bin_coverage_median(cov, bin_size = 10)
    ## 3.2 Process coverage
    min.limit <- min( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "start"] )
    max.limit <- max( gtf[gtf$type == "transcript" & gtf$transcript_id %in% tx.annot.df$transcript_id, "end"] )
    cov <- process.cov.interactionDTU(cov, 
                                    condition.levels = c(condition_1, condition_2), 
                                    condition.col = "POP_Condition", # for interaction DTU we average coverage per POP_Condition
                                    min.limit,
                                    max.limit, 
                                    type = "mean")
    ## 3.3 Coverage Track
    groups.i <- setdiff(colnames(cov), c("chr", "start", "end"))
    ordered.groups <- c("AFB_NS", "AFB_IAV", "AFB_COV", "EUB_NS", "EUB_IAV", "EUB_COV") 
    groups.i <- ordered.groups[ordered.groups %in% groups.i]
    # You must order the coverage. The legend will correspond with the order of the columns in cov
    cov <- cov[, c("chr", "start", "end", groups.i)]
    cov.track <-  coverage.track.fun(cov,
                                     groups = groups.i, 
                                     groups.color = pop.condition.color.vec[ groups.i ], 
                                     track.name = "Mean coverage")
    
    # 4. Merge
    track.l <- c(tx.track.l, cov.track )
    
    # 5. Equivalence class Track
    ## 5.1 get EC coordinates (from the entire EC set)
    ec.coord.l <- EC.coordinate.generator( gtf, gene_id, ecs = all.ecs )
    ec.names <- names(ec.coord.l)
    ec.coord.l <- lapply( ec.names, function (ec) {
      ec.coord <- ec.coord.l[[ec]]
      if(length(ec.coord) > 0){  ec.coord$feature_name <- ec } # add feature name to GRanges
      return(ec.coord)
    })
    names(ec.coord.l) <- ec.names
    ec.coord.gtf <- Reduce(c, ec.coord.l)
    
    ## 5.2 ECs tested in DTU will be coloured and those not tested will be in grey or some other color outside of the pallette
    tested.ecs <- ec.annot.df$transcript_name
    
    ## 5.3 highlight tracks must be overlayed on top of an existing track list (tx.track and coverage track in our case)
    track.l <- EC.highlight.fun( tracklist = track.l, ec.coord.gtf )
    ## 5.4 Color tracks
    track.l <- edit.EC.highlight.track.fill( track.l, ec.annot.df )
    return(track.l)
    
    # close tryCatch
    }, error = function(e) e) 
  })
  # add names
  gene_name.vec <- tx2gene[ match(dtu.df[["gene_id"]], tx2gene$gene_id), "gene_name"]
  names(track.l) <- gene_name.vec
  # Remove NULL entries (failed iterations)
  track.l <- Filter(Negate(is.null), track.l)
  return(track.l)
}