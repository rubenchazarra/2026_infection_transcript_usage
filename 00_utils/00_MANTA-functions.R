#!/usr/bin/env Rscript

##  0. Packages ------------------------------------------

suppressPackageStartupMessages(require(dplyr))
# For MANTA
suppressPackageStartupMessages(require(manta))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(reshape2))
# For MANTARRAYA: MANTA edit to run with paired data
suppressPackageStartupMessages(require(car))
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(CompQuadForm))

##  1. Processing  ------------------------------------------

group_by_Condition_Aquino <- function(counts.l){
  # Group list of Condition count dataframes ready for pairwise comparison: "COV-NS", "IAV-NS", "COV-IAV"
  groups <- list("IAV.NS" = c("IAV", "NS"), "COV.NS" = c("COV", "NS"), "IAV.COV" = c("IAV", "COV"))
  
  if(class(counts.l[[1]]) == "list"){ # if object inside list is list, iterate. This is for when we have data split by POP
    lapply(counts.l, group_by_Condition)
  }else{
    counts.l <- lapply(groups, function(x) {
      l <- list(counts.l[[ x[1] ]], counts.l[[ x[2] ]])
      setNames(l, nm = x)
    })
    names(counts.l) <- names(groups)
    return(counts.l)
  }
}

add.gene_id.to.transcript_name <- function(counts){
  ## This function adds gene id as column to a df with transcript names as rownames. Adding gene_id and transcript_id
  
  tx2gene <- read.table("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/01_Quantification/00_IndexData/01_GENCODEv42/gencode.v42-IAV-tx2gene-with-ENS.tsv", header = T, sep = "\t", na.strings = "")
  tx2gene <- tx2gene[ -c( grep("PAR", tx2gene$gene_id) ), ]
  ## checks
  print("table(rownames(counts) %in% tx2gene$transcript_name) : ...")
  table(rownames(counts) %in% tx2gene$transcript_name)
  counts <- merge(counts, tx2gene[, c("gene_id", "transcript_id", "transcript_name")], by.x=0, by.y = "transcript_name", all.x = T, all.y = F)
  rownames(counts) <- counts[["Row.names"]]
  counts <- counts[, -c(1)]
  counts$transcript_name <- NULL
  # Order columns
  ord.cols <- c("gene_id", "transcript_id", names(counts)[grep("HMN", names(counts))])
  return(counts[, ord.cols])
}

annotate.equivalence.classess <- function(ec.vec, tx2gene, n_cores){
  # print("Adding gene_name, gene_id, transcript_names to vector of ECs ...")
  ec.vec <- gsub("-", "_", ec.vec )
  # Annotate ECs with a single transcript_id (matching works here)
  annot.df = cbind( tx2gene[ match(ec.vec, tx2gene$transcript_id), c("gene_id", "gene_name", "transcript_name") ], 
                    data.frame("ec" = ec.vec ))
  # Annotate ECs representing > 1 transcript_id
  ec.idx <- grep(",", annot.df[["ec"]])
  if(length(ec.idx) > 0){
    ec <- annot.df[ec.idx, "ec"]
    ec.l <- parallel::mclapply(ec, function(x){
      idx <- match(unlist(strsplit(x, ",")), tx2gene$transcript_id)
      df = data.frame( "gene_name" = paste0( unique(tx2gene[idx, "gene_name"]), collapse = ","),
                       "gene_id" = paste0( unique(tx2gene[idx, "gene_id"] ), collapse = ","), 
                       "transcript_name" = paste0( sort(unique(tx2gene[idx, "transcript_name"] )),
                                                   collapse = ","))
    }, mc.cores = 12)
    ec.df <- plyr::ldply(ec.l)
    annot.df[ec.idx, "gene_name"] <- ec.df$gene_name
    annot.df[ec.idx, "gene_id"] <- ec.df$gene_id
    annot.df[ec.idx, "transcript_name"] <- ec.df$transcript_name
    annot.df[["ec"]] <- ec.vec
  }
  rownames(annot.df) <- annot.df[["ec"]]
  return(annot.df)
}

sQTLseekeR2.rel_abundances <- function(counts, min.transcript.exp=1, min.gene.exp=5, min.prop=0.4, min.dispersion= 0.01) {
  ## Compute relative abundance of transcripts per gene in counts
  print("Computing EC relative abundance ...")
  annot.cols <- c("gene_id", "gene_name", "transcript_name", "ec")
  annot.ec.df <- counts[, annot.cols]
  # Edit names for sQTLSeeker2 requirements
  names(counts)[ match(c("gene_id", "ec"), names(counts)) ] <- c("geneId", "trId")
  ## Remove unneeded columns
  counts <- counts[, -c ( match(c("gene_name", "transcript_name"), names(counts)) )]
  
  ## 5. Calculate relative expression (sQTLseekeR2::prepare.trans.exp())
  rel_counts <- sQTLseekeR2::prepare.trans.exp(te.df = counts, 
                                               min.transcript.exp = min.transcript.exp, 
                                               min.gene.exp = min.gene.exp, 
                                               min.prop = min.prop, 
                                               min.dispersion = min.dispersion, # 0.1 # main parameter in filtering
                                               verbose = F)
  # NOTE: NAs are introduced in counts samples with lower gene expression than min.gene.exp
  print(paste0("Number of genes passing filtering: ", length(unique(rel_counts[["geneId"]])), " ..."))
  # Edit colnames back to original 
  names(rel_counts)[ match(c("geneId", "trId"), names(rel_counts)) ] <- c("gene_id", "ec")
  # Add gene_name and transcript_name back 
  rel_counts <- cbind( rel_counts, annot.ec.df [ match(rel_counts[["ec"]], annot.ec.df[["ec"]]),  c( "gene_name", "transcript_name")] )
  rel_counts <- rel_counts %>% relocate(gene_name, .after = "gene_id" )
  rel_counts <- rel_counts %>% relocate(transcript_name, .after = "ec" )
  rownames(rel_counts) <- rel_counts[["ec"]]
  return(rel_counts)
}

sQTLseekeR2.rel_abundances.wrapper <- function(counts, tx2gene, meta_data, n_cores, min.transcript.exp=1, min.gene.exp=5, min.prop=0.8, min.dispersion= 0.01, meta_col ){
  if( class(counts) == "list"){  lapply(counts, function(counts) sQTLseekeR2.rel_abundances.wrapper(counts, tx2gene, meta_data, n_cores, min.transcript.exp, min.gene.exp, min.prop, min.dispersion, meta_col)) 
  } else {
    
    # 1. Annotate ECs (add gene_id, gene_name)
    annot.ec.df <- annotate.equivalence.classess(ec.vec = rownames(counts), tx2gene, n_cores = 112 )
    counts <- cbind( annot.ec.df, counts )
    
    # 2. Compute relative abundances
    rel_counts <- sQTLseekeR2.rel_abundances( counts, min.transcript.exp, min.gene.exp, min.prop, min.dispersion)  # if fails, try running LANG=C.UTF-8 # before run
    annot.cols <- c( "gene_id", "gene_name", "transcript_name", "ec")
    annot.ec.df <- rel_counts[ , annot.cols ]
    
    # 2. Compute MD (effect size) with rel_counts
    donor_cond.vec <-  colnames(counts)[ -match(annot.cols, colnames(counts)) ]
    condition.vec <- meta_data[ match( donor_cond.vec, meta_data[["Donor.ID_Condition"]] ), meta_col ]
    md.df <- md.trans.wrapper( rel_counts, annot.cols, condition.vec, n_cores )
      
    # 3. Calculate highest mean EC per gene per Condition or Population
    donor_cond.vec.l <- split( donor_cond.vec, condition.vec )
    condition.l <- lapply(donor_cond.vec.l, function(x) rel_counts[ , match(x, colnames(rel_counts))] )
    ec.switch.df <- mean.ec.per.condition ( condition.l, annot.ec.df[, c("gene_id", "gene_name", "transcript_name")] )
    ec.switch.df[["condition_1"]] <- names(condition.l)[1]
    ec.switch.df[["condition_2"]]  <- names(condition.l)[2]
    
    # 4. Merge outputs
    annot.ec.df <- merge( annot.ec.df, ec.switch.df, by = c("gene_id", "gene_name"), all = T)
    annot.ec.df <- merge( annot.ec.df, md.df, by = "gene_id", all = T)
    rel_counts = merge(annot.ec.df, rel_counts, by = annot.cols, all = T) 
    rownames(rel_counts) <- rel_counts[["ec"]] # rownames are shuffled in merge()
    return(rel_counts)
  }
}

shuffle.donor.labels <- function(rel_df, index){
  # Shufle relative_abundance donor labels based on precomputed shuffled label file
  print(paste0("*** Running permutation [ ", index, " ] ***"))
  # 2. Shuffle meta.data labels
  label.file <- "/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/01_Randolph2021/05_Differential-Analyses/01_DEA/00_labels-permutations/01_label-permutation-90-donors.txt"
  label.df <- read.table( label.file, header = F)
  new.labels <- as.character( label.df[ index, ])
  new.labels <- new.labels[ new.labels %in% meta_data[["infection_ID"]]]
  ## re-label donors
  id.cols <- c("gene_id", "gene_name", "tcc", "transcript_name")
  idx <- match(setdiff(colnames(rel_df), id.cols), colnames(rel_df))
  colnames(rel_df)[idx] <- new.labels # this does not work
  # reorder counts labels as in metadata
  rel_df <- counts_edits(rel_df, meta_data)
  return(rel_df)
}

md.trans <- function(sr.o, groups.o) {
  # sQTLseekeR2 function to calculate absolute maximum difference in mean adjusted EC relative abundance between conditions (MD value)  (https://github.com/guigolab/sQTLseekeR2/blob/master/R/md.trans.R)
  mTrans <- apply(sr.o, 2, function(sr.r) tapply(sr.r, groups.o, mean, na.rm = TRUE))
  lr <- nrow(mTrans)
  ind1 <- rep(1:(lr-1), (lr-1):1)
  ind2 <- NULL
  for(ii in 2:lr){
    ind2 <- c(ind2, ii:lr)
  }
  MDtrans <- apply(mTrans, 2, function(r) diff(rbind(r[ind1], r[ind2])))
  if(!is.matrix(MDtrans)){
    tr.names <- names(MDtrans)
    MDtrans <- matrix(MDtrans, 1)
    colnames(MDtrans) <- tr.names
  }
  gpMD <- apply(MDtrans, 1, function(e) max(abs(e)))
  gpMD.max <- which.max(gpMD)
  tr.first <- names(which.max(abs(MDtrans[gpMD.max, ])))
  tr.second <- names(which.max(-sign(MDtrans[gpMD.max, tr.first]) * MDtrans[gpMD.max, ]))
  df <- data.frame("MD" =  max(gpMD, na.rm = TRUE), "tr.first" = tr.first, "tr.second" = tr.second)
  return(df)
}

md.trans.wrapper <- function( counts, annot.cols, condition.vec, n_cores ){
  # Calculate maximum difference (MD) in relative abundance between conditions. Effect size estimate implemented in sQTLseekeR2::md.trans()
  md.l <- parallel::mclapply( unique(counts$gene_id), function( gene_id ){
    i <- counts[counts$gene_id %in% gene_id, ]
    rownames(i) <- i[["transcript_name"]] # This is to have transcript_name as id for MD df
    i <- t( i[ ,  - match( annot.cols, colnames(i))])
    md.trans( sr.o = i, groups.o = condition.vec )
  }, mc.cores = n_cores)
  
  md.df <- data.table::rbindlist(md.l)
  md.df[["gene_id"]] <- unique(counts$gene_id)
  return(md.df)
}

mean.ec.per.condition <- function(condition.l, annot.df){
  print("Computing most abundant EC per condition ...")
  
  mean.l <- lapply(condition.l, function(x) cbind(annot.df, data.frame(mean_ec_exp = rowMeans(x, na.rm = T))) )
  mean.l <- lapply(mean.l, function(x) {
    x %>% group_by(gene_id, gene_name) %>% slice_max(mean_ec_exp)
  })
  ec.switch.df <- merge(mean.l[[1]], mean.l[[2]], by = c("gene_id", "gene_name"), all = T)
  names(ec.switch.df) <- c ( "gene_id", "gene_name", "max.tx_name.1",  "mean_ec_exp.1", "max.tx_name.2", "mean_ec_exp.2" )
  ec.switch.df$is.EC.switch <- ec.switch.df$max.tx_name.1 !=  ec.switch.df$max.tx_name.2
  return(ec.switch.df)
}


##  2. MANTA (independent data DTU) ------------------------------------------
run_manta <- function(gene_id, relative_abundances, meta_data, covariates, annot_cols, transform="none", ss.type="III"){
  # Run MANTA for DTU analysis
  ## Format data
  input_data <- t(relative_abundances[ , which(! names(relative_abundances) %in%  annot_cols )]) # remove annotation columns
  formula = as.formula(paste("input_data ~ ", paste(c(covariates), collapse = " + ")))
  ## Run MANTA
  manta_obj <- manta(
    formula = formula,
    data = meta_data,
    transform = transform, 
    type = ss.type) # previously using "I", @ Diego Garrido (method author) suggests using SS type = "II" / "III" for DTU analysis
  ## Extract MANTA output 
  manta_df = as.data.frame(manta_obj$aov.tab[ covariates , ])
  ## Rename output df columns
  colnames(manta_df) <- c("Df","Sum Sq", "Mean Sq", "F value", "R2", "p_value")
  manta_df[["feature_ids"]] <- paste0(colnames(input_data), collapse=";")
  manta_df[["n_features"]] <- ncol(input_data)
  manta_df[["covariate"]] <- rownames(manta_df)
  ## add gene metadata 
  annot_vec <- apply(relative_abundances[, annot_cols] , 2, function(x) paste0(unique(x), collapse=";"))
  annot_df = t(data.frame(annot_vec))
  rownames(annot_df) <- NULL
  # add
  manta_df <- cbind(annot_df, manta_df)
  # add N of samples the test is performed with and N of samples excluded
  manta_df[["samples_tested"]] <- table(apply(input_data, 1, function(x) all(!is.na(x))))[["TRUE"]]
  manta_df[["samples_excluded"]] <- nrow(input_data) - manta_df[["samples_tested"]]
  return(manta_df)
}

p_value_correction <- function(df, p_value.col = "p_value"){
  # Benjamini-Hochberng multiple test correction of p_values
  df[["fdr"]] <- p.adjust(df[[p_value.col]], method = "BH")
  rownames(df) <- df[["gene_id"]]
  return(df)
}

print_manta_stats <- function(manta_results.l){
  # Print stats of MANTA execution for each of the covariates in the input model
  print.l <- lapply(names(manta_results.l), function(df_name){
    df = manta_results.l[[df_name]]
    data.frame("TRUE" = ifelse("TRUE" %in% names(table(df$fdr < 0.05)), yes = table(df$fdr < 0.05)["TRUE"], no = 0), 
               "FALSE" = table(df$fdr < 0.05)["FALSE"], row.names = df_name )
  })
  print.df <- Reduce(rbind, print.l)
  print("==================")
  print(print.df)
  print("==================")
}

manta_wrapper <- function(rel_abundances, meta_data, covariates, meta_col = "Donor.ID_Condition") {
  # Iterativelly run MANTA for differential transcript abundance (DTU) analysis considering input covariates present in meta_data
  
  if(class(rel_abundances) == "list") {lapply(rel_abundances, function(rel_abundances) manta_wrapper(rel_abundances, meta_data, covariates ))
    
  }else{ 
    print("Running MANTA ...")
    
    # 1. Subset and order metadata to match rel_abundances samples 
    meta_col <- "Donor.ID_Condition" #  ifelse( any ( colnames(rel_abundances) %in%  meta_data[["Sample.ID"]] ), yes = "Sample.ID", no = "Donor.ID_Condition" )
    sample_cols <- intersect(colnames(rel_abundances),  meta_data[[meta_col]] )
    annot_cols <-  setdiff(colnames(rel_abundances),  sample_cols )
    meta_data <- meta_data[ match(sample_cols, meta_data[[meta_col]]), ]
    progress.vec <- unique(rel_abundances$gene_id) [ seq(0, length(unique(rel_abundances$gene_id)), by = 100) ] 
    
    # 2. If all values in one metadata column are the same, remove covariate
    covariates <- covariates[ sapply(meta_data[, covariates], function(x) length(unique(x)) > 1) ]
    model <- as.formula(paste("dtu ~ ", paste(c(covariates), collapse = " + ")))
    print( model )
    
    # 3. Run DTU 
    manta_results <- lapply(unique(rel_abundances$gene_id),  function(gene_id) { 
      ## Print gene index
      if(gene_id %in% progress.vec ) print(match(gene_id, unique(rel_abundances$gene_id)))
      ## Run MANTA
      tryCatch({ # to deal with potential errors
        run_manta(gene_id, 
                  relative_abundances = rel_abundances[rel_abundances$gene_id %in% gene_id, ], 
                  meta_data, 
                  covariates, 
                  annot_cols,
                  transform="none", 
                  ss.type = "III")
      }, error = function(e) {
        message("Error in run_manta for gene: ", gene_id)
        message(" -> ", e$message)
        return(NULL)
      })  
    }) #,  mc.cores = n_cores, mc.preschedule = FALSE) # using parallel::mclapply leads to job freezing and never ending
    
    # 3. Process output
    ## 3.1 Remove errors # @RCG 2025-03-20 This was sorted a while ago
    keep.dfs <- unlist(lapply(manta_results, is.data.frame))
    print(table(keep.dfs))
    manta_results <- manta_results[keep.dfs]
    manta_results <- as.data.frame(data.table::rbindlist(manta_results))
    ## 3.2 Split by covariate
    manta_results.l <- split( manta_results, manta_results[["covariate"]])
    
    # 4. Multiple test correction
    manta_results.l <- lapply(manta_results.l, function(x) p_value_correction(df = x, p_value.col = "p_value"))
    # 5. Print stats
    print_manta_stats(manta_results.l)
    
    return(manta_results.l)
  }
}

manta_permutations_wrapper <- function(rel_abundances,  meta_data, covariates, n_permutations) { 
  # Run MANTA permutations 
  ## Read shuffled labels
    shuffled.labels <- readRDS("../../../../00_Misc/01_labels-permutations/01_populationDTU-shuffle-labels-x100.rds")
    annot.cols <- c("gene_id", "gene_name", "transcript_name", "ec" )
    
    N_it <- n_permutations # 25  # 2024-07-01 running now 10 iterations, I will scale up later
    iteration_l <- lapply(c( 1: N_it ), function(it){
      print(paste0("MANTA permutation N: ", it, " .................."))
      
      # 1. Rename labels
      rel_abundances <- mapply(function(rel_df, label_df){
        sample_cols <- intersect( colnames(rel_df),  meta_data[["Donor.ID_Condition"]] )
        annot_cols <-  setdiff(colnames(rel_df),  sample_cols )
        sample.df <- rel_df[, sample_cols ]
        annot.df <- rel_df[, annot_cols ]
        annot.df <- annot.df[,  c("gene_id", "gene_name", "transcript_name", "ec" )] # remove gene metadata which is not interesting because we are going to recompute DTU with shuffled labels
        # Add shuffled names
        shuffled.ids <- unlist(label_df[it, , drop = T])
        itsc.names <- intersect( shuffled.ids, colnames(rel_df)) # this is for the case where not all samples are present e.g at the cell type level.
        sample.df <- setNames(sample.df, itsc.names) # rename
        cbind(annot.df, sample.df)
      }, rel_abundances, shuffled.labels, SIMPLIFY = F)
      
      # 2.Run MANTA
      manta_wrapper(rel_abundances, meta_data, covariates, meta_col = "Donor.ID_Condition")
    })
    ## add names
    names( iteration_l ) <- paste0("iteration_", c(1 : N_it) )
    return(iteration_l)
}

##  3. MANTARRAYA (paired data DTU) ------------------------------------------

## This is a set of function termed MANTARRAYA developed by Diego Garrido, Ferran Reverter and Miquel Calvo to conduct DTU with paired data
## The input data is the substraction of vectors of relative abundances between 2 conditions, say rel_abundances (infected) - rel_abundances (baseline)
##
remove.most.corr.tx <- function(tx_data){
  # Given a matrix where columns are tx and rows are samples, remove the most correlated tx of the matrix with another tx of the matrix
  corr.mat = cor(tx_data)
  corr.mat[lower.tri(corr.mat, diag = T)] <- 0 # remove lower triangle and diagonal
  # remove transcripts with correlations > 0.75 btw them
  max.corr.tx = names(which.max(apply(abs(corr.mat), 2, max)))
  tx_data = tx_data[, grep(max.corr.tx, colnames(tx_data), invert = TRUE), drop = F] # remove most corr tx
  return(tx_data)
}

mantarraya_more2tx = function(fit, dc){
  # This function is a wrapper to run mantarraya which thought to run DTU on paired data samples. For instance relative abundance differences between 2 conditions of the same sample. E.g changes upon infection. 
  # The interpretation of the results is slightly different to standard manta.
  # In mantarraya a model such as relative abundance (infected - non infected) = Intercept + age, here the intercept tells if there are significant differences upon infection, and beta age is the interaction of infection_status * age
  # Developed by Diego Garrido Martin (2024-02)
  # @RCG 2024-05-13 MANTA for genes with > 2 Tx 
  UU = car::Anova(fit, type = "III") # sums of squares
  SS = unlist(lapply(UU$SSP, function(x){ sum(diag(x)) })) 
  SS.e = sum(diag(UU$SSPE))
  
  # df = c(1, fit$assign[-1])   # degrees of freedom # change suggested by @Diego
  df = as.numeric(table(fit$assign)) # degrees of freedom
  
  df.e = fit$df.residual 
  
  f.tilde = SS/SS.e * df.e/df # test statistics 
  
  lambda = eigen(cov(fit$residuals), only.values = T)$values
  lambda <- lambda[lambda/sum(lambda) > 10^-3] 
  
  # Compute p-values & print summary stats
  acc = 1e-14
  tests = c("Intercept", names(dc))
  # Create a data frame with the results
  out.df = data.frame(test = tests, SS = SS, pseudo.F = f.tilde, df = df, row.names = tests)
  p.value.vec <- c()
  for ( i in 1:length(tests)){   
    p.val =  farebrother(q = SS[i], lambda, h = rep(df[i], length(lambda)), eps = acc)$Qq
    p.value.vec <- c(p.value.vec, p.val)
  }
  out.df[["p.value"]]<- p.value.vec
  return(out.df)
}

mantarraya_2tx <- function(fit, dc){
  # @RCG 2024-05-13 MANTARRAYA (relative abunadnce test of paired data) for the case of genes with 2 Tx
  UU = car::Anova(fit, type = "III")  
  f.tilde = UU$`F value`   # estadísticos (intercept, factores y un NA al final)
  df = UU$Df        # grados de libertad (intercept, factores y el residuo al final)
  SS = UU[["Sum Sq"]] # sum of squares

  tests = c("Intercept", names(dc), "Residuals") 
  out.df = data.frame(test = tests, SS = SS, pseudo.F = f.tilde, df = df, row.names = tests)
  p.value.vec <- c()
  for ( i in 1:(length(df)-1) ) {
    p.val <- pchisq( f.tilde[i]*df[i] , df[i] , lower.tail = F)  # p-value for the i-th term
    p.value.vec <- c(p.value.vec, p.val)
  }
  p.value.vec <- c( p.value.vec, NA ) # add NA for residuals term
  out.df[["p.value"]]<- p.value.vec
  out.df[["feature_ids"]] <- paste0(colnames(fit$coefficients), collapse=";")
  return(out.df)
}

# Run MANTARRAYA 
run_mantarraya_edit <- function( input_data, covariates, meta_data){
  # Run mantarraya distinguishing set ups depending on the N of Tx
  model <- if( all(covariates == "")){ as.formula(paste("input_data ~ 1")) }else{ as.formula(paste("input_data ~ ", paste(c(covariates), collapse = " + "))) }
  fit <- lm( model, data = meta_data )
  dc <- attr(terms(fit),"dataClasses")[-1]
  
  n_tx <- ncol(input_data)
  if( n_tx == 1 ){ # Mantarraya for genes with 2 Tx
    mantarraya_2tx(fit, dc)
  } else if( n_tx >= 2) { # Mantarraya for genes with >= 2 Tx
    mantarraya_more2tx(fit, dc)
  }
}

mantarraya_loop <- function( input_data, covariates, meta_data ) {
  # Run MANTARRAYA in an iterative manner, due to close to total anticorrelation of the transcript ratio vectors, this can lead to errors, if this happens we remove one 1 Tx and model 1 instead of 2
  obj_list <- tryCatch({
    manta_obj <- run_mantarraya_edit ( input_data, covariates, meta_data ) # RUN MANTARRAYA
    failed <- FALSE # update failed status if reaching here (this is, not failing above)
    list(manta_obj, failed, input_data) # out
    }, error = function(e) { # print("Mantarraya failed, removing last TCC/Tx ...")  # stop printing, if not std.out files get huge
    manta_obj <- NULL
    failed <- TRUE
    input_data <- input_data[, -c( ncol(input_data) ), drop = F ] # remove last Tx from input_data # TODO: remove less expressed Tx
    list(manta_obj, failed, input_data) # return updated failed status, and input_data
  })
  # Update values
  failed = obj_list[[2]] ; input_data = obj_list[[3]]
  # if MANTARRAYA failed, re run with updated values
  if( failed ){ 
    obj_list <- mantarraya_loop(input_data, covariates, meta_data) 
  }else{
    mr_out = obj_list[[1]]
    mr_out[["feature_ids"]] <- paste0(colnames(input_data), collapse=";")
    mr_out[["n_features"]] <- length( colnames(input_data))
    return(mr_out)
  }
}

mantarraya_wrapper <- function(gene_id, relative_abundances, meta_data, covariates, annot_cols){
  # This function is a wrapper to run MANTARRAYA: an edit of the MANTA DTU test to work on paired data samples.
  # Given inputing paired to MANTA breaks the requirement of independence of samples, we input differences in relative abundances between conditions
  # For instance here we model the difference in relative abundance between donor_A_infected - donor_A_Baseline 
  # The interpretation of the results is slightly different to standard MANTA 
  # In MANTARRAYA a model such as relative abundance (infected - non infected) = Intercept + age, here the intercept tells if there are significant differences upon infection, and beta age is the interaction of infection_status * Age
  
  # 1) Relative abundance matrix of gene_id
  input_data <- t ( relative_abundances[, which(!names(relative_abundances) %in% annot_cols)] )
  
  # 2) Run MANTARRAYA
  mr_out <- mantarraya_loop( input_data, covariates, meta_data )
  
  # 3) Add Additional gene metadata
  # gene_meta <- relative_abundances[ relative_abundances$gene_id %in% gene_id, annot_cols]
  # gene_meta <- unique(gene_meta[ , - match (c("transcript_name", "ec"), colnames(gene_meta))]) # transcript_name and ec columns contain 1 value per row and we are collapsing to one per gene. We remove them
  # rownames(gene_meta) <- NULL 
  # mr_out <- merge(gene_meta, mr_out)
  
  # Add transcripts excluded in mantarraya_loop()
  tx.excluded <- setdiff(colnames(input_data), unlist(strsplit(mr_out[["feature_ids"]], ";")))
  n.tx.excluded <- length(tx.excluded)
  mr_out[["n_features.tested"]] <- mr_out[["n_features"]]
  mr_out[["n_features.excluded"]] <- n.tx.excluded
  mr_out[["n_features"]] <- mr_out[["n_features"]] + mr_out[["n_features.excluded"]]
  mr_out[["feature_ids.excluded"]] <- paste0(tx.excluded, collapse=";")
  ## add gene metadata
  gene_meta <- relative_abundances[ relative_abundances$gene_id %in% gene_id, annot_cols]
  annot_vec <- apply(relative_abundances[, annot_cols] , 2, function(x) paste0(unique(x), collapse=";"))
  annot_df = t(data.frame(annot_vec))
  rownames(annot_df) <- NULL
  # add 
  mr_out <- cbind(annot_df, mr_out)
  # add N of samples the test is performed with and N of samples excluded
  mr_out[["samples_tested"]] <- table(apply(input_data, 1, function(x) all(!is.na(x))))[["TRUE"]]
  mr_out[["samples_excluded"]] <- nrow(input_data) - mr_out[["samples_tested"]]
  return(mr_out)
}

mantarraya_meta_wrapper <- function(rel_abundances, meta_data, covariates ){ 

  if(class(rel_abundances) == "list") { lapply(rel_abundances, function(rel_abundances) mantarraya_meta_wrapper(rel_abundances, meta_data, covariates))
    
  }else{ 
    condition_1 <- unique(rel_abundances$condition_1)
    condition_2 <- unique(rel_abundances$condition_2)
    print(paste0("Running MANTARRAYA for Condition : ", condition_1, " - ", condition_2, " ..."))
    # 1. Subset and order metadata to match rel_abundances samples 
    ## 2025-08-06 Subset metadata to the first condition (e.g for IAV_NS keep IAV). This selection is done to select cell type abundances # 2025-08-30 For some reason condition_1 is now "AFB" ??? Bug appeared while running the filter benchmark
    # meta_data <- meta_data[meta_data$Condition %in% condition_1, ]
    meta_col <- "Donor.ID" # difference in relative abudances matrix contains Donor.ID in colnames in Aquino et al., 2023 dataset
    sample_cols <- intersect(colnames(rel_abundances),  meta_data[[meta_col]] )
    annot_cols <-  setdiff(colnames(rel_abundances),  sample_cols )
    meta_data <- meta_data[ match(sample_cols, meta_data[[meta_col]]), ]
    progress.vec <- unique(rel_abundances$gene_id) [ seq(0, length(unique(rel_abundances$gene_id)), by = 100) ] 
    
    # 2. If all values in one metadata column are the same, remove covariate
    covariates <- covariates[ sapply(meta_data[, covariates], function(x) length(unique(x)) > 1) ]
    model <- as.formula(paste("dtu ~ ", paste(c(covariates), collapse = " + ")))
    print( model )
    
    # 3. Run MANTARRAYA DTU 
    mr_results <- lapply( unique(rel_abundances$gene_id), function(gene_id) { 
      # Progress bar
      if(gene_id %in% progress.vec ) print(match(gene_id, unique(rel_abundances$gene_id)))
      tryCatch({ # to deal with potential errors
        ## Run MANTARRAYA
        mantarraya_wrapper ( gene_id, 
                             relative_abundances = rel_abundances [ rel_abundances$gene_id %in% gene_id, ], 
                             meta_data, 
                             covariates, 
                             annot_cols )
      }, error = function(e) {
        message("Error in run_manta for gene: ", gene_id)
        message(" -> ", e$message)
        return(NULL)
      })
    })# , mc.cores = 56, mc.preschedule = FALSE) # 2025-05-06 running parallel::mclapply inside parallel::mclapply will lead to error
    
    # 3. Process output
    ## 3.1. Remove errors
    keep.dfs <- unlist(lapply(mr_results, is.data.frame))
    mr_results <- mr_results[keep.dfs]
    mr_results <- as.data.frame(data.table::rbindlist(mr_results))
    ## 3.2. Split by covariate
    mr_results.l <- split( mr_results, mr_results[["test"]])
    
    # 4. Multiple test correction
    mr_results.l <- lapply(mr_results.l, function(x) p_value_correction(df = x, p_value.col = "p.value"))
    # 5. Print stats
    print_manta_stats(mr_results.l)
    return(mr_results.l)
  }
}

mantarraya_permutations_wrapper <- function(rel_abundances, meta_data, covariates, n_permutations) { 
  
  # 1. Processing 
  ## 1.1 Remove unwanted gene metadata columns
  annot_cols <- c("gene_id", "gene_name", "transcript_name", "ec")
  
  # 2. Shuffled labels
  shuffled.labels <- readRDS("/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/03_Differential/00_Misc/01_labels-permutations/02_03_infectionDTU-shuffle-infection-labels-x100.rds")
  shuffled.labels <- shuffled.labels[c("COV.NS", "IAV.NS")] # 2025-10-01 Interested in COV-Baseline and IAV-Baseline. Not in IAV COV comparison. 
  # 3. Iterate across Conditions
  condition.iteration.l <- lapply( names(shuffled.labels), function(cond){
    print(paste0("Running iteration for Condition: ", cond, " .................."))
    
    rel_abundances <- rel_abundances[[cond]]
    # 3.1. Keep donors present in both conditions
    meta_col <- "Donor.ID_Condition"
    sample_cols <- intersect(colnames(rel_abundances),  meta_data[[meta_col]] )
    annot_cols <-  setdiff(colnames(rel_abundances),  sample_cols )
    annot_df <- rel_abundances[, annot_cols]
    rel_abundances <- rel_abundances[, sample_cols]
    meta_data <- meta_data[ match(sample_cols, meta_data[[meta_col]]), ]
    
    ## Subset rel_abundances to Donors samples present in both Conditions (for the posterior subtraction)
    meta.i <- meta_data[, c("Donor.ID", "Donor.ID_Condition", "Condition")]
    count.df <- meta.i %>% group_by(Donor.ID) %>% summarise(n_conditions = n_distinct(Condition))
    samples.select <- meta.i[meta.i$Donor.ID %in% count.df[count.df$n_conditions > 1, "Donor.ID", drop = T], "Donor.ID_Condition"]
    rel_abundances <- rel_abundances[, samples.select]
    
    # 3.4 Iterate across permutations
    cond.labels <- shuffled.labels[[ cond ]]
    N_it <- n_permutations # 10 
    iteration.l <- lapply(c(1:N_it), function(i){
      
      print(paste0("Iteration N: ", i, " .................."))
      # 1. Split rel_abundance matrix into 2 according to shuffled labels
      labels.1 <-  sort( intersect(colnames(rel_abundances), as.character( cond.labels[[1]][i, ])) ) # sort to substract donors in a paired manner
      labels.2 <-  sort( intersect(colnames(rel_abundances), as.character(cond.labels[[2]][i, ]))  )
      ## Check donor ordering is the same 
      table(gsub("_NS|_COV|_IAV", "", labels.1) == gsub("_NS|_COV|_IAV", "", labels.2))
      # split counts by shuffled labels
      mat.l <- list("matrix_1" = rel_abundances[ , labels.1 ], "matrix_2" = rel_abundances[ , labels.2 ] )
      mat.l <- lapply(mat.l, function(mat) setNames(mat, gsub("_NS|_COV|_IAV", "", colnames(mat))))
      # 2. Substract matrices (MANTARRAYA set up)
      counts = mat.l[[1]] - mat.l[[2]] 
      counts = cbind( annot_df, counts) 
      rownames(counts) <- counts[["ec"]]
      
      # 3. Run MANTA
      mantarraya_meta_wrapper( rel_abundances = counts, meta_data, covariates )
    })
    names(iteration.l) <- paste0("permutation_", c(1:N_it)) # add names
    return(iteration.l)
  })
  names(condition.iteration.l) <- names(shuffled.labels) # add names
  return(condition.iteration.l)
}