#### FUNCTIONS for general Data Analysis of the Single-cell Alternative Splicing project #####

read.RDS.files <- function(path, pattern=".rds"){
  # Read RDS files from a directory into a list and name them
  file.vec <- list.files(path, pattern = pattern, full.names = T)
  file.l <- lapply(file.vec, readRDS)
  setNames(object = file.l, nm = gsub("-.*", "", basename(file.vec)))
}

read.txt.files <- function(path, pattern=".txt"){
  # Read RDS files from a directory into a list and name them
  file.vec <- list.files(path, pattern = pattern, full.names = T)
  file.l <- lapply(file.vec, function(x) read.table(x, sep = "\t", header = T))
  setNames(object = file.l, nm = gsub("-.*", "", basename(file.vec)))
}

add.gene_name.2.gene_id <- function(df, tx2gene, gene_id.col, gene_name.col){
  print("Adding gene_name to gene_id")
  # Add gene name to gene_id
  if(class(df) == "list"){ lapply(df, function(df) add.gene_name.2.gene_id( df, tx2gene, gene_id.col, gene_name.col ))
  }else{
    df[[gene_name.col]] <- tx2gene[ match(df[[gene_id.col]], tx2gene$gene_id), "gene_name"]
    return(df)
  }
}

add.tx_name.2.tcc_id <- function(df, tx2gene, tcc_id.col, tx_name.col ){ # FIISH THIS FUNCTION!!!!S
  print("Converting TCCs from transcript_id to transcript_name")
  # Convert TCC names to tx names, sort by N of tx included per TCC
  if(class(df) == "list"){ lapply(df, function(df) add.tx_name.2.tcc_id( df, tx2gene, tcc_id.col, tx_name.col ))
    }else{
      tx_name.l <- lapply(strsplit(df[[tcc_id.col]], ";"), function(tcc_id.vec){
        tx_name.l <- lapply(strsplit(tcc_id.vec, ","), function(tx_id.vec){
          tx_name_vec <- tx2gene[match(tx_id.vec, tx2gene$transcript_id) , "transcript_name"]
          paste(tx_name_vec, collapse = ",")
          })
        paste(unlist(tx_name.l), collapse = ";")
        })
      df[[tx_name.col]] <- unlist(tx_name.l)
        
  return(df)
    }
}


filter.fdr <- function(df, fdr.col = "fdr", fdr.thres =0.05){
  ## Filter a dataframe by fdr.col{
  if(class(df) == "list"){ lapply(df, function(x) filter.fdr(x, fdr.col))
  } else {
    df = df[df[[fdr.col]] < fdr.thres ,  ]
    print(paste0("N DE feat: ", nrow(df)))
    return(df)
  }
}

filt.FC.and.fdr <- function(l, FC.col, FC.thres, fdr.col, fdr.thres){
  if(class(l) == "list") lapply(l, function(l) filt.FC.and.fdr( l, FC.col, FC.thres, fdr.col, fdr.thres ) )
  else{
    # filter by abs log2FC and FDR
    l <- l[ abs(l [, FC.col] ) > FC.thres & l [, fdr.col] < fdr.thres, ]
    print(paste0("N DE feat: ", nrow(l)))
    return(l)
  }
}

label.DEG <- function(de, FC.col, FC.thres, fdr.col, fdr.thres){
  if(class(de) == "list") lapply(de, function(de) label.DEG(de, FC.col, FC.thres, fdr.col, fdr.thres ) )
  else{
    # filter by abs log2FC and FDR
    de <- de %>%
      mutate(
        is.DGE = .data[[fdr.col]] < fdr.thres & abs(.data[[FC.col]]) > FC.thres,
        DGE.direction = case_when(
          .data[[fdr.col]] < fdr.thres & .data[[FC.col]] > FC.thres  ~ "Up",
          .data[[fdr.col]] < fdr.thres & .data[[FC.col]] < -FC.thres ~ "Down",
          TRUE ~ "Not DGE"
          ))
    print(paste0("N DE feat: ", nrow(de[de$is.DGE, ])))
    return(de)
  }
}

filt.DEG <- function(de, is.dge.col){
  if(class(de) == "list") lapply(de, function(de) filt.DEG(de, is.dge.col ) )
  else{
    de <- de[ de[[is.dge.col]], ] # filt DE
    return(de)
  }
}


select.df.from.nested.list <- function(l, df.name){
   # Select particular item from nested list of dataframes
   if( class(l[[1]]) == "list" ){ lapply(l, function(l) select.df.from.nested.list(l , df.name))
   }else{
       l[[ grep(df.name, names(l)) ]]
     }
 }


n.dge.up.down <- function(df, FC.col="betas"){
  if(class(df) == "list") lapply(df, function(df) n.dge.up.down(df, FC.col))
  else{
    data.frame("n_dge" = c( nrow(df),
                            nrow(df[df[, grep(FC.col, names(df))] > 0, ]),
                            nrow(df[df[, grep(FC.col, names(df))] < 0, ])),
               "Class" = c("All", "Up", "Down")
    )
  }
}


test.overlap <- function(m1, m2, fdr.col = "fdr", gene.col = "gene_id", FC = FALSE, FC.col = NULL, FC.thres=0.5){
  # Get significance of overlap of 2 models with "gene_id" column
  if(FC){
    sig_1 = m1 [ m1[[fdr.col]] < 0.05 & abs(m1[[FC.col]]) > FC.thres , gene.col ]
    sig_2 = m2 [ m2[[fdr.col]] < 0.05 & abs(m2[[FC.col]]) > FC.thres , gene.col ]
  }else{
    sig_1 = m1 [ m1[[fdr.col]] < 0.05 , gene.col ]
    sig_2 = m2 [ m2[[fdr.col]] < 0.05 , gene.col ]
  }
  
  all_1 = m1[[gene.col]]
  all_2 = m2[[gene.col]]
  
  df <- data.frame("shared" = length( intersect( sig_1, sig_2 )), 
                   "model.1.specific" = length( setdiff( sig_1, sig_2 )), # model specific but SIG
                   "model.2.specific" = length( setdiff( sig_2, sig_1 )) )
  # contingency table
  cont.tab <- contingency.table.fun(sig_1, sig_2, all_1, all_2)
  ## add contingency table values
  df[["S12"]] <- cont.tab[1, 1] ; df[["S1_NS2"]] <- cont.tab[2, 1] ; df[["S2_NS1"]] <- cont.tab[1, 2] ; df[["NS12"]] <- cont.tab[2, 2]
  
  # Fisher exact test
  fisher_result <- fisher.test(cont.tab)
  df[["p_value"]] <- fisher_result$p.value
  df[["odds.ratio"]] <- fisher_result$estimate[["odds ratio"]]
  df[["conf.int.lower"]] <- fisher_result$conf.int[1]
  df[["conf.int.upper"]] <- fisher_result$conf.int[2]
  
  return(df)
}

contingency.table.fun <- function( s1, s2, all_1, all_2 ){
  
  ns1 <- setdiff(all_1, s1)
  ns2 <- setdiff(all_2, s2)
  
  s12 <- intersect(s1, s2)
  s1_ns2 <- intersect(s1, ns2)
  s2_ns1 <-  intersect(s2, ns1)
  ns12 <- intersect(all_1, all_2)
  
  contingency_table <- matrix(c(length(s12), # significant in 1 & 2
                                length(s1_ns2), # NS in 1, S in 2
                                length(s2_ns1), # NS in 2, S in 1 
                                length(ns12) - length(s12)), ncol = 2, byrow = TRUE) # not significant in both
  return(contingency_table)
}


test.overlap.pro <- function(m1, m2, 
                             fdr.col.1, fdr.thres.1, gene.col.1, FC.col.1 = NA, FC.thres.1 = 0, 
                             fdr.col.2, fdr.thres.2, gene.col.2, FC.col.2 = NA, FC.thres.2 = 0){
  if(class(m1)[1] == "list" & class(m2)[1] == "list"){ mapply( function(m1, m2) test.overlap.pro(m1, m2, 
                                                                         fdr.col.1, fdr.thres.1, gene.col.1, FC.col.1, FC.thres.1, 
                                                                         fdr.col.2, fdr.thres.2, gene.col.2, FC.col.2, FC.thres.2), m1, m2, SIMPLIFY = F)
  }else{
  # Get significance of overlap of 2 models with "gene_id" column
  if(!is.na(FC.col.1) &  !is.na(FC.col.2) ){ # 1. Two models with logFC & FDR
    print("A")
    sig_1 = m1 [ m1[[fdr.col.1]] < fdr.thres.1 & abs(m1[[FC.col.1]]) > FC.thres.1 , gene.col.1 , drop = T]
    sig_2 = m2 [ m2[[fdr.col.2]] < fdr.thres.2 & abs(m2[[FC.col.2]]) > FC.thres.2 , gene.col.2 , drop = T]
  }else if(is.na(FC.col.1) & !is.na(FC.col.2)){ # 2. Model1 without FC, model2 with FC (e.g comparing DTU vs DGE)
    print("B")
    sig_1 = m1 [ m1[[fdr.col.1]] < fdr.thres.1 , gene.col.1 , drop = T]
    sig_2 = m2 [ m2[[fdr.col.2]] < fdr.thres.2 & abs(m2[[FC.col.2]]) > FC.thres.2 , gene.col.2 , drop = T]
  }else if(!is.na(FC.col.1) & is.na(FC.col.2)){ # 3. Model1 with FC, model2 without FC (e.g comparing DGE vs DTU )
    print("C")
    sig_1 = m1 [ m1[[fdr.col.1]] < fdr.thres.1 & abs(m1[[FC.col.1]]) > FC.thres.1 , gene.col.1 , drop = T]
    sig_2 = m2 [ m2[[fdr.col.2]] < fdr.thres.2 , gene.col.2 , drop = T]
  }else if(is.na(FC.col.1) & is.na(FC.col.2)){ # 4. Two models with FDR only (e.g DTU vs DTU)
    print("D")
    sig_1 = m1 [ m1[[fdr.col.1]] < fdr.thres.1 , gene.col.1, drop = T]
    sig_2 = m2 [ m2[[fdr.col.2]] < fdr.thres.2 , gene.col.2, drop = T]
  }
  
  all_1 = m1[[gene.col.1]]
  all_2 = m2[[gene.col.2]]
  
  df <- data.frame("shared" = length( intersect( sig_1, sig_2 )), 
                   "model.1.specific" = length( setdiff( sig_1, sig_2 )), # model specific but SIG
                   "model.2.specific" = length( setdiff( sig_2, sig_1 )) )
  # build contingency table
  cont.tab <- contingency.table.fun(sig_1, sig_2, all_1, all_2)
  
  ## add contingency table values
  df[["S12"]] <- cont.tab[1, 1] ; df[["S1_NS2"]] <- cont.tab[2, 1] ; df[["S2_NS1"]] <- cont.tab[1, 2] ; df[["NS12"]] <- cont.tab[2, 2]
  
  # Fisher exact test
  fisher_result <- fisher.test(cont.tab)
  df[["p_value"]] <- fisher_result$p.value
  df[["odds.ratio"]] <- fisher_result$estimate[["odds ratio"]]
  df[["conf.int.lower"]] <- fisher_result$conf.int[1]
  df[["conf.int.upper"]] <- fisher_result$conf.int[2]
  
  return(df)
  }
}


test.fisher <- function(sig_1, sig_2, all_1, all_2){
  # contingency table
  cont.tab <- contingency.table.fun(sig_1, sig_2, all_1, all_2)
  ## add contingency table values
  df <- data.frame("T1T2" = NA)
  df[["T1T2"]] <- cont.tab[1, 1] ; df[["T1F2"]] <- cont.tab[2, 1] ; df[["T2F1"]] <- cont.tab[1, 2] ; df[["F1F2"]] <- cont.tab[2, 2]
  
  # Fisher exact test
  fisher_result <- fisher.test(cont.tab)
  df[["p_value"]] <- fisher_result$p.value
  df[["odds.ratio"]] <- fisher_result$estimate[["odds ratio"]]
  df[["conf.int.lower"]] <- fisher_result$conf.int[1]
  df[["conf.int.upper"]] <- fisher_result$conf.int[2]
  df <- add.p_value.label(df, p_value.col ="p_value") # add p value label
  return(df)
}

add.p_value.label <- function(df, p_value.col ="p_value"){
  print( "Add labels for pvalue in ggplot ..." )
  if(class(df)[1] == "list"){ lapply(df, function(df) add.p_value.label(df, p_value.col))
  }else{
    df$p_value_label <- ""
    df$p_value_label [ df[[p_value.col]] > 0.05 ] <-'ns'
    df$p_value_label [ df[[p_value.col]] < 0.05 ] <-'*'
    df$p_value_label [ df[[p_value.col]] < 1e-02 ] <- '**' # 0.01
    df$p_value_label [ df[[p_value.col]] < 1e-03 ] <- '***'
    df$p_value_label [ df[[p_value.col]] < 1e-04 ] <- '****'
    df$p_value_label [ df[[p_value.col]] < 1e-05 ] <- '*****'
    return(df)
  }
}

select.first.N  <- function(df, sort.col, N){
  # Select first N elements from df after sorting bby sort.col
  if(class(df) == "list"){ lapply(df, function(df) select.first.N(df, sort.col, N)) # iterate if list
  }else{
    df <- df[order( df[[sort.col]], decreasing = F ), ]
    if( nrow(df) > N ) df <- df[ c(1:N) , ]
    return(df)
  }
}

get.nrow <- function(df){
  if(class(df) == "list"){ lapply(df, get.nrow )
  }else{
    data.frame("n" = nrow(df))
  }
}

prepare.ComplexUpset <- function(vec.l){
  print("Prepraing input list for ComplexUpset plot ...")
  # remove empty
  vec.l <- vec.l[ unlist(lapply(vec.l, length)) > 0 ]
  unique_elements <- unique(unlist(vec.l))
  result_df <- data.frame(matrix(FALSE, nrow = length(unique_elements), 
                                 ncol = length(vec.l)))
  rownames(result_df) <- unique_elements
  colnames(result_df) <- names(vec.l)
  # Step 3: Populate the dataframe with 1's where elements are present
  for (vec_name in names(vec.l)) {
    result_df[vec.l[[vec_name]], vec_name] <- TRUE
  }
  return(result_df)
}

# prepare.ComplexUpset.2 <- function(vec.l){
#   print("Prepraing input list for ComplexUpset plot ...")
#   # remove empty
#   vec.l <- vec.l[ unlist(lapply(vec.l, length)) > 0 ]
#   unique_elements <- unique(unlist(vec.l))
#   result_df <- data.frame(matrix(FALSE, 
#                                  ncol = length(unique_elements), 
#                                  nrow = length(vec.l)))
#   colnames(result_df) <- unique_elements
#   rownames(result_df) <- names(vec.l)
#   # Step 3: Populate the dataframe with 1's where elements are present
#   for (vec_name in names(vec.l)) {
#     result_df[vec_name, vec.l[[vec_name]]] <- TRUE
#   }
#   return(result_df)
# }
