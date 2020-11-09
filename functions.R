##########
#
# Nacer et al.
# Pan-cancer application of a lung-adenocarcinoma-derived 
# gene-expression-based prognostic predictor
#
# Functions used in the main script
#
##########

library(tidyverse)
library(org.Hs.eg.db)


# Loading data  -----

loadRData <- function(fileName){
  # loads an RData file and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Changing row names of matrices  -----

# from reporter id to...
gex_matrix_reporter_to_symbol <- function(gex_matrix, gene_table, reporter_col, symbol_col, dup_removal_fun) {
  # uses gene tables that come with the matrix
  # changes row names from reporter ids to gene symbols
  # deals with duplicates according to the dup_removal_fun function
  # very slow code, needs to be improved...
  gene_symbols <- gex_matrix[,1:2]
  gene_symbols <- rownames_to_column(data.frame(gene_symbols), var = reporter_col)
  gene_symbols <- left_join(gene_symbols, gene_table[,c(reporter_col,symbol_col)])
  gex_matrix <- as.data.frame(gex_matrix)
  gex_matrix$geneSymbol <- gene_symbols[[symbol_col]] # add symbol column to gex_matrix
  gex_matrix <- gex_matrix %>% group_by(geneSymbol) %>% summarise_all(dup_removal_fun) # remove all duplicates
  gex_matrix <- gex_matrix[complete.cases(gex_matrix[[symbol_col]]),] # remove all NAs from the symbol column
  gex_matrix <- column_to_rownames(data.frame(gex_matrix), var = symbol_col) # rename the rows of the matrix
  gex_matrix <- as.matrix(gex_matrix)
  return(gex_matrix)
}


gex_matrix_reporter_to_symbol2 <- function(gex_matrix, gene_table, reporter_col, symbol_col, dup_removal_fun) {
  # uses gene tables that come with the matrix
  # changes row names from reporter ids to gene symbols
  # deals with duplicates according to the dup_removal_fun function
  # faster than previous code
  gene_symbols <- gex_matrix[,1:2]
  gene_symbols <- rownames_to_column(data.frame(gene_symbols), var = reporter_col) # have a column of gene symbols
  gene_symbols <- left_join(gene_symbols, gene_table[,c(reporter_col,symbol_col)]) # same order as reporters
  gex_matrix <- as.data.frame(gex_matrix)
  chunks <- c(seq(0, length(colnames(gex_matrix)), 100), length(colnames(gex_matrix))) # divide matrix in blocks of 100 columns + remainder
  for (i in head(seq_along(chunks), -1)) {
    begin_slice <- chunks[i]+1 # e.g. 1, 101, 201...
    end_slice <- chunks[i+1] # e.g. 100, 200, 300...
    matrix_chunk <- gex_matrix[,begin_slice:end_slice] # slice the matrix
    matrix_chunk$geneSymbol <- gene_symbols[[symbol_col]] # add symbol column to gex_matrix
    matrix_chunk <- matrix_chunk %>% group_by(geneSymbol) %>% summarise_all(dup_removal_fun) # remove all duplicates
    matrix_chunk <- matrix_chunk[complete.cases(matrix_chunk[[symbol_col]]),] # remove all rows with NA in the symbol column
    ifelse(exists("gex_matrix_reduced"), # save all reduced into a new dataframe
           gex_matrix_reduced <- left_join(gex_matrix_reduced, matrix_chunk, by=symbol_col),
           gex_matrix_reduced <- matrix_chunk)
  }
  gex_matrix_reduced <- column_to_rownames(gex_matrix_reduced, var = symbol_col) # rename the rows of the matrix
  gex_matrix_reduced <- as.matrix(gex_matrix_reduced)
  return(gex_matrix_reduced)
}


names_reporter_to_symbol <- function(gex_matrix, gene_table, reporter_col, symbol_col) {
  # uses gene tables that come with the matrix
  # gets a list of gene symbols in the same order as the reporters in the matrix
  reporters <- gex_matrix[,1:2]
  reporters <- rownames_to_column(data.frame(reporters), var = reporter_col)
  reporters <- left_join(reporters, gene_table[,c(reporter_col,symbol_col)])
  gene_symbols <- reporters[[symbol_col]]
  return(gene_symbols)
}

gex_matrix_reporter_to_entrez <- function(gex_matrix, gene_table, reporter_col, entrez_col, dup_removal_fun, add_e) {
  # uses gene tables that come with the matrix
  # changes row names from reporter ids to entrez ids (possibly with an 'e' in front)
  # deals with duplicates according to the dup_removal_fun function
  entrez_ids <- gex_matrix[,1:2]
  entrez_ids <- rownames_to_column(data.frame(entrez_ids), var = reporter_col)
  entrez_ids <- left_join(entrez_ids, gene_table[,c(reporter_col,entrez_col)])
  gex_matrix <- as.data.frame(gex_matrix)
  gex_matrix$entrezId <- entrez_ids[[entrez_col]] # add symbol column to gex_matrix
  gex_matrix <- gex_matrix %>% group_by(entrezId) %>% summarise_all(dup_removal_fun) # remove all duplicates
  gex_matrix <- gex_matrix[complete.cases(gex_matrix$entrezId),] # remove all NAs from the symbol column
  if (add_e == 'yes_add_e') {
    gex_matrix <- gex_matrix %>% mutate(entrezId = paste0("e", entrezId))
  }
  gex_matrix <- column_to_rownames(data.frame(gex_matrix), var = "entrezId") # rename the rows of the matrix
  gex_matrix <- as.matrix(gex_matrix)
  return(gex_matrix)
}


# from ensembl id to...
gex_matrix_ensembl_to_symbol <- function(gex_matrix, gene_table, ensembl_col, symbol_col) {
  # uses gene tables that come with the matrix
  # changes row names from ensembl ids to gene symbols
  ensembl_ids <- gex_matrix[,1:2]
  ensembl_ids <- rownames_to_column(data.frame(ensembl_ids), var = ensembl_col)
  ensembl_ids <- left_join(ensembl_ids, gene_table[,c(ensembl_col,symbol_col)])
  gene_symbols <- ensembl_ids[[symbol_col]]
  row.names(gex_matrix) <- gene_symbols
  return(gex_matrix)
}


gex_matrix_ensembl_to_symbol_Hs <- function(gex_matrix) {
  # uses org.Hs.eg.db, for when gene information table is not available
  # changes row names from ensembl ids to gene symbols
  ensembl_ids <- row.names(gex_matrix)
  gene_symbols <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column='SYMBOL', keytype='ENSEMBL')
  gene_symbols <- unname(gene_symbols)
  row.names(gex_matrix) <- gene_symbols
  return(gex_matrix)
}

gex_matrix_ensembl_to_entrez <- function(gex_matrix, gene_table, ensembl_col, entrez_col, add_e) {
  # uses gene tables that come with the matrix
  # changes row names from ensembl ids to gene symbols
  ensembl_ids <- gex_matrix[,1:2]
  ensembl_ids <- rownames_to_column(data.frame(ensembl_ids), var = ensembl_col)
  ensembl_ids <- left_join(ensembl_ids, gene_table[,c(ensembl_col,entrez_col)])
  entrez_ids <- ensembl_ids[[entrez_col]]
  if (add_e == 'yes_add_e') {
    entrez_ids <- paste0("e",entrez_ids)
  }
  row.names(gex_matrix) <- entrez_ids
  return(gex_matrix)
}

gex_matrix_ensembl_to_entrez_Hs <- function(gex_matrix, add_e) {
  # uses org.Hs.eg.db, for when gene information table is not available
  # changes row names from ensembl ids to entrez ids (possibly with an 'e' in front)
  ensembl_ids <- row.names(gex_matrix)
  entrez_ids <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column='ENTREZID', keytype='ENSEMBL')
  entrez_ids <- unname(entrez_ids)
  if (add_e == 'yes_add_e') {
    entrez_ids <- paste0("e",entrez_ids)
  }
  row.names(gex_matrix) <- entrez_ids
  return(gex_matrix)
}


# from gene symbol to...
gex_matrix_symbol_to_entrez_Hs <- function(gex_matrix, add_e) {
  # uses org.Hs.eg.db, for when gene information table is not available
  # changes row names from reporter ids to entrez ids (possibly with an 'e' in front)
  gene_symbols <- row.names(gex_matrix)
  entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_symbols, column='ENTREZID', keytype='SYMBOL')
  entrez_ids <- unname(entrez_ids)
  if (add_e == 'yes_add_e') {
    entrez_ids <- paste0("e",entrez_ids)
  }
  row.names(gex_matrix) <- entrez_ids
  return(gex_matrix)
}

names_symbol_to_entrez_Hs <- function(gex_matrix) {
  # uses org.Hs.eg.db, for when gene information table is not available
  # gets a list of entrez ids in the same order as the gene symbols in the matrix
  gene_symbols <- row.names(gex_matrix)
  entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_symbols, column='ENTREZID', keytype='SYMBOL')
  entrez_ids <- unname(entrez_ids)
  return(entrez_ids)
}


# For CLAMS  -----

# get_samples_class_from_predictor_result
get_class_from_aims_result <- function(aims_result) {
  # get sample's classification from AIMS result
  sample.class <- aims_result$cl[,]
  sample.class <- data.frame(as.list(sample.class))
  sample.class <- t(sample.class)
  sample.class <- rownames_to_column(data.frame(sample.class), var = "sample.id")
  return(sample.class)
}


get_probability_from_aims_result <- function(aims_result) {
  # get sample's classification probability from AIMS result
  sample.prob <- aims_result$prob[,]
  sample.prob <- data.frame(as.list(sample.prob))
  sample.prob <- t(sample.prob)
  sample.prob <- rownames_to_column(data.frame(sample.prob), var = "sample.id")
  return(sample.prob)
}


count_duplicated_genes <- function(list_with_duplicates) {
  # count duplicated genes
  uniques <- unique(list_with_duplicates)
  dup_genes <- c()
  for (i in 1:length(uniques)) {
    gene.count <- str_count(string=list_with_duplicates, pattern=regex(paste0("^", uniques[i], "$"))) %>% sum()
    dup_genes <- c(dup_genes, gene.count)
  }
  return(dup_genes)
}


get_all_genes_in_rules <- function(all.pairs){
  # get all genes present in the rules (keeps duplicates)
  genes <- c()
  for (cp in strsplit(all.pairs,"<")){
    genes <- c(genes,cp)
  }
  return(genes)
}


get_all_unique_genes_in_rules <- function(all.pairs){
  # get all unique genes from the rules of an SSP
  genes <- get_all_genes_in_rules(all.pairs)
  return(unique(genes))
}


get_rules_per_gene <- function(all.pairs) {
  # count in how many rules a gene appears in
  rules.genes <- get_all_genes_in_rules(all.pairs)
  unique.genes <- get_all_unique_genes_in_rules(all.pairs)
  genes_in_rules <- c()
  for (i in 1:length(unique.genes)) {
    rule.count <- str_count(string=rules.genes, pattern=regex(paste0(unique.genes[i], "$"))) %>% sum()
    genes_in_rules <- c(genes_in_rules, rule.count)
  }
  return(genes_in_rules)
}


fix_symbols_for_clams <- function(gex_matrix) {
  # input: a gene expression matrix with gene symbols as row names
  # and samples as column names
  # outputs the same matrix with the gene symbols as CLAMS uses them
  gene_names <- row.names(gex_matrix)
  fixed_names <- case_when(gene_names == 'VEGFD' ~ 'FIGF',
                           gene_names == 'H2AZ1' ~ 'H2AFZ',
                           gene_names == 'CAVIN2' ~ 'SDPR',
                           TRUE ~ gene_names)
  row.names(gex_matrix) <- fixed_names
  return(gex_matrix)
}



# Gene signature  -----

get_signature_value_for_matrix <- function(gex_matrix, gene_list, dup_removal_fun) {
  # returns a value for a gene signature (gene_list)
  # for all samples (columns) from a gene expression matrix
  # dealing with duplicates with the function passed as dup_removal_fun
  #
  # create empty frame for result
  sample.id <- colnames(gex_matrix)
  results <- as.data.frame(sample.id)
  gex_matrix <- rownames_to_column(gex_matrix, var='geneSymbol')
  geneSymbol <- gex_matrix["geneSymbol"]
  # get value for every column that is not the geneSymbol one
  for (i in 2:length(colnames(gex_matrix))) {
    sample_column <- gex_matrix[,i,drop=FALSE]
    # create a subset with gene symbols and sample
    current.sample <- bind_cols(data.frame(geneSymbol), data.frame(sample_column))
    # remove duplicate rows by getting max value for each symbol (if dup_removal_fun = max)
    nodup <- current.sample %>% group_by(geneSymbol) %>% summarise_all(dup_removal_fun)
    # order by value
    ordered <- nodup %>% arrange(.[[2]])
    # extract position for list of genes (module), sum them and add to results table
    module.rows <- which(ordered$geneSymbol %in% gene_list)
    module.value <- sum(as.integer(module.rows))
    results$value[results$sample.id == colnames(sample_column)[1]] <- module.value
    rm(sample_column)
    if (i %% 50 == 0) { print(i) }
  }
  return(results)
}

# For overall survival analysis
convert_time <- function (patient_table) {
  # mutate patient OS info from original months/years to months and years
  patient_table$OS.time.type <- as.character(patient_table$OS.time.type)
  patient_table$OS.time <- as.double(patient_table$OS.time)
  patient_table <- patient_table %>%
    mutate(OS.time.years = case_when(OS.time.type == 'days' ~ (OS.time / 365.2425),
                                     OS.time.type == 'months' ~ (OS.time / 12),
                                     OS.time.type == 'years' ~ OS.time,
                                     TRUE ~ OS.time),
           OS.time.months = case_when(OS.time.type == 'days' ~ (OS.time / 30.436875),
                                      OS.time.type == 'months' ~ OS.time,
                                      OS.time.type == 'years' ~ (OS.time * 12),
                                      TRUE ~ OS.time))
  return(patient_table)
}

censor_data <- function (patient_table, time_to_censor, unit_to_censor) {
  # censor OS column to a specific time point
  # get new names and what columns to use for censoring
  time_column_to_use <- case_when(unit_to_censor == "y" ~ "OS.time.years",
                                  unit_to_censor == "m" ~ "OS.time.months")
  censored_time_column <- case_when(unit_to_censor == "y" ~ paste0("OS.time.years.", time_to_censor, unit_to_censor),
                                    unit_to_censor == "m" ~ paste0("OS.time.months.", time_to_censor, unit_to_censor))
  censored_event_column <- paste0("OS.event.", time_to_censor, unit_to_censor)
  
  # censor
  patient_table <- patient_table %>%
    # if OS number > time_to_censor, set time to cap, else keep original number
    mutate(censored_time =
             ifelse(get(time_column_to_use) > as.numeric(time_to_censor),
                    as.numeric(time_to_censor), get(time_column_to_use)),
           # if OS number > time_to_censor, set event to 0, else keep original event
           censored_event = 
             ifelse(get(time_column_to_use) > as.numeric(time_to_censor),
                    0, OS.event))
  names(patient_table)[names(patient_table) == "censored_time"] <- censored_time_column
  names(patient_table)[names(patient_table) == "censored_event"] <- censored_event_column  
  return(patient_table)
}
