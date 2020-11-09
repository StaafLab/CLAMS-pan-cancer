##########
#
# Nacer et al.
# Pan-cancer application of a lung-adenocarcinoma-derived 
# gene-expression-based prognostic predictor
#
# Reference code
# for reproducing analysis after collecting data
# and generating figures
#
##########

# load packages, config file and functions
library(tidyverse)
library(yaml)
library(Biobase)
library(CLAMSstrict)
library(cowplot)

config <- yaml.load_file("config_datasets.yml")

source(config$functions_path)

# Classification with CLAMS  ------
# Separate analysis for each data set

# GOBO ________________________________________________________________
# input data
setwd(config$data_directory_path)
gex_matrix <- loadRData(config$gobo$gex_matrix$file)
gene_table <- loadRData(config$gobo$gene_table$file)
reporter_col <- config$gobo$gene_table$reporter_col
symbol_col <- config$gobo$gene_table$symbol_col

# convert rows to gene symbols
gex_matrix_gobo <- gex_matrix_reporter_to_symbol2(gex_matrix, gene_table, reporter_col, symbol_col, 'max')
rm(gene_table)
rm(gex_matrix)

# run CLAMS and extract classification results
gex_matrix_gobo_reduced <- filter_gex_for_ssp(gex_matrix_gobo, CLAMSmodel)
clams_result_gobo <- applyCLAMS(gex_matrix_gobo_reduced, CLAMSmodel)
clams.class <- get_class_from_aims_result(clams_result_gobo)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'GOBO'

# save results to file
setwd(config$result_dir_path)
write.csv(clams.class, "gobo_clams.csv", row.names=FALSE)


# SCAN-B ________________________________________________________________
# input data
setwd(config$data_directory_path)
# matrix already with gene symbols
gex_matrix_scanb <- loadRData(config$scanb$gex_matrix$file)

# run CLAMS and extract classification results
ssp_genes <- get_all_unique_genes_in_rules(CLAMSmodel$all.pairs)
gex_matrix_scanb_reduced <- gex_matrix_scanb[ssp_genes,]
clams_result_scanb <- applyCLAMS(gex_matrix_scanb_reduced, CLAMSmodel)
clams.class <- get_class_from_aims_result(clams_result_scanb)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'SCAN-B'

# save results to file
setwd(config$result_dir_path)
write.csv(clams.class, "scanb_clams.csv", row.names=FALSE)


# TCGA ________________________________________________________________
# input data
setwd(config$data_directory_path)
gex_matrix <- exprs(readRDS(config$tcga$rds$file))

# convert rows to gene symbols
gex_matrix_tcga <- gex_matrix_ensembl_to_symbol_Hs(gex_matrix)
rm(gex_matrix)
gex_matrix_tcga <- fix_symbols_for_clams(gex_matrix_tcga) # fix gene symbols if needed

# run CLAMS and extract classification results
ssp_genes <- get_all_unique_genes_in_rules(CLAMSmodel$all.pairs)
gex_matrix_tcga_reduced <- gex_matrix_tcga[ssp_genes,]
clams_result_tcga <- applyCLAMS(gex_matrix_tcga_reduced, CLAMSmodel)
clams.class <- get_class_from_aims_result(clams_result_tcga)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'TCGA'

# save results to file
setwd(config$result_dir_path)
write.csv(clams.class, "tcga_clams.csv", row.names=FALSE)


# Tabchy GSE20271 ________________________________________________________________
# input data
setwd(config$data_directory_path)
gex_matrix <- loadRData(config$gse20271$gex_matrix$file)
gene_table <- loadRData(config$gse20271$gene_table$file)
reporter_col <- config$gse20271$gene_table$reporter_col
symbol_col <- config$gse20271$gene_table$symbol_col

# convert rows to gene symbols
gex_matrix_gse20271 <- gex_matrix_reporter_to_symbol2(gex_matrix, gene_table, reporter_col, symbol_col, 'max')
rm(gene_table)
rm(gex_matrix)

# run CLAMS and extract classification results
clams_result_gse20271 <- run_clams(gex_matrix_gse20271)
clams.class <- get_class_from_aims_result(clams_result_gse20271)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'GSE20271'

# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "gse20271_clams.csv", row.names=FALSE)


# Iwamoto GSE22093 ________________________________________________________________
# input data
setwd(config$data_directory_path)
gex_matrix <- loadRData(config$gse22093$gex_matrix$file)
gene_table <- loadRData(config$gse22093$gene_table$file)
reporter_col <- config$gse22093$gene_table$reporter_col
symbol_col <- config$gse22093$gene_table$symbol_col

# convert rows to gene symbols
gex_matrix_gse22093 <- gex_matrix_reporter_to_symbol2(gex_matrix, gene_table, reporter_col, symbol_col, 'max')
rm(gene_table)
rm(gex_matrix)

# run CLAMS and extract classification results
clams_result_gse22093 <- run_clams(gex_matrix_gse22093)
clams.class <- get_class_from_aims_result(clams_result_gse22093)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'GSE22093'

# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "gse22093_clams.csv", row.names=FALSE)


# Hess ________________________________________________________________
# input data
setwd(config$data_directory_path)

gex_matrix <- read_delim(config$hess$gex_matrix$file, delim = '\t')
symbol_col <- config$hess$gex_matrix$symbol_col

# remove columns that are not samples
# and convert the column of gene symbols to row names
gex_matrix$reporterId <- NULL
gex_matrix$entrezId <- NULL
gex_matrix_hess <- gex_matrix[complete.cases(gex_matrix[,symbol_col]),] # remove rows that have no gene symbol
gex_matrix_hess <- column_to_rownames(gex_matrix_hess, symbol_col)
gex_matrix_hess <- as.matrix(gex_matrix_hess)
rm(gex_matrix)

# run CLAMS and extract classification results
clams_result_hess <- run_clams(gex_matrix_hess)
clams.class <- get_class_from_aims_result(clams_result_hess)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'Hess'

# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "hess_clams.csv", row.names=FALSE)


# Melanoma
# Riaz BMS038 ________________________________________________________________
setwd(config$data_directory_path)
gex_matrix_bms038 <- read.delim(config$bms038$gex_matrix$file,as.is=T)
gex_matrix_bms038 <- gex_matrix_bms038[-which(is.na(gex_matrix_bms038$HUGO)),]
table(duplicated(gex_matrix_bms038$HUGO)) # no duplicates
row.names(gex_matrix_bms038) <- NULL
gex_matrix_bms038 <- column_to_rownames(gex_matrix_bms038, var="HUGO")
gex_matrix_bms038 <- as.matrix(gex_matrix_bms038)

# run CLAMS and extract classification results
clams_result_bms038 <- run_clams(gex_matrix_bms038)
clams.class <- get_class_from_aims_result(clams_result_bms038)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'BMS038'
# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "bms038_clams.csv", row.names=FALSE)


# Liu ________________________________________________________________
setwd(config$data_directory_path)

gex_matrix_schadendorf <- read_delim(config$schadendorf$gex_matrix$file, delim='\t')
gex_matrix_schadendorf <- column_to_rownames(gex_matrix_schadendorf, var="X1") # patients as rows, gene symbols as columns
gex_matrix_schadendorf <- t(gex_matrix_schadendorf)

# run CLAMS and extract classification results
clams_result_schadendorf <- run_clams(gex_matrix_schadendorf)
clams.class <- get_class_from_aims_result(clams_result_schadendorf)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'Schadendorf'
# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "schadendorf_clams.csv", row.names=FALSE)


# Gide ________________________________________________________________
setwd(config$data_directory_path)

gex_matrix_gide <- read_delim(config$gide$gex_matrix$file, delim='\t')
table(duplicated(gex_matrix_gide$X1)) # no duplicates
gex_matrix_gide <- column_to_rownames(gex_matrix_gide, var='X1')
gex_matrix_gide <- as.matrix(gex_matrix_gide)

# run CLAMS and extract classification results
clams_result_gide <- run_clams(gex_matrix_gide)
clams.class <- get_class_from_aims_result(clams_result_gide)
clams.class <- clams.class %>% dplyr::rename(clams.class = sample.class)
clams.class["dataset"] <- 'Gide'
# save results to file
setwd(config$result_dir_path)

write.csv(clams.class, "gide_clams.csv", row.names=FALSE)



# Proliferation values  -----
# Find genes in common in all data sets
genes_gobo <- row.names(gex_matrix_gobo)
genes_scanb <- row.names(gex_matrix_scanb)
genes_tcga <- row.names(gex_matrix_tcga)
genes_gse20271 <- row.names(gex_matrix_gse20271)
genes_gse22093 <- row.names(gex_matrix_gse22093)
genes_hess <- row.names(gex_matrix_hess)
genes_bms038 <- row.names(gex_matrix_bms038)
genes_schadendorf <- row.names(gex_matrix_schadendorf)
genes_gide <- row.names(gex_matrix_gide)

common_genes_9 <- Reduce(intersect, list(genes_gobo, genes_scanb, genes_tcga,
                                  genes_gse20271, genes_gse22093, genes_hess,
                                  genes_bms038, genes_schadendorf, genes_gide))

# input proliferation signature module from Karlsson et al (93 genes)
mod_prolif_karlsson <- read_csv("modules/Karlsson/Karlsson.csv")
mod_prolif_karlsson <- mod_prolif_karlsson$geneSymbol

prolif_karlsson_in_common_9 <- mod_prolif_karlsson[mod_prolif_karlsson %in% common_genes_9]

common_genes_list <- common_genes_9
prolif_karlsson_in_common <- prolif_karlsson_in_common_9

# GOBO ________________________________________________________________
# filter matrix to common genes
gex_matrix_gobo <- gex_matrix_gobo[row.names(gex_matrix_gobo) %in% common_genes_list,]
gex_matrix_gobo <- as.data.frame(gex_matrix_gobo)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_gobo, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "gobo_karlsson_9.csv", row.names=FALSE)


# SCAN-B ________________________________________________________________
gex_matrix_scanb <- gex_matrix_scanb[row.names(gex_matrix_scanb) %in% common_genes_list,]
gex_matrix_scanb <- as.data.frame(gex_matrix_scanb)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_scanb, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "scanb_karlsson_9.csv", row.names=FALSE)


## TCGA ________________________________________________________________
gex_matrix_tcga <- gex_matrix_tcga[row.names(gex_matrix_tcga) %in% common_genes_list,]
gex_matrix_tcga <- as.data.frame(gex_matrix_tcga)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_tcga, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "tcga_karlsson_9.csv", row.names=FALSE)


# GSE20271 ________________________________________________________________
gex_matrix_gse20271 <- gex_matrix_gse20271[row.names(gex_matrix_gse20271) %in% common_genes_list,]
gex_matrix_gse20271 <- as.data.frame(gex_matrix_gse20271)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_gse20271, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "gse20271_karlsson_9.csv", row.names=FALSE)


# GSE22093 ________________________________________________________________
gex_matrix_gse22093 <- gex_matrix_gse22093[row.names(gex_matrix_gse22093) %in% common_genes_list,]
gex_matrix_gse22093 <- as.data.frame(gex_matrix_gse22093)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_gse22093, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "gse22093_karlsson_9.csv", row.names=FALSE)


# Hess ________________________________________________________________
gex_matrix_hess <- gex_matrix_hess[row.names(gex_matrix_hess) %in% common_genes_list,]
gex_matrix_hess <- as.data.frame(gex_matrix_hess)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_hess, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "hess_karlsson_9.csv", row.names=FALSE)


# Riaz ________________________________________________________________
gex_matrix_bms038 <- gex_matrix_bms038[row.names(gex_matrix_bms038) %in% common_genes_list,]
gex_matrix_bms038 <- as.data.frame(gex_matrix_bms038)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_bms038, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "bms038_karlsson_9.csv", row.names=FALSE) # change name accordingly


# Liu ________________________________________________________________
gex_matrix_schadendorf <- gex_matrix_schadendorf[row.names(gex_matrix_schadendorf) %in% common_genes_list,]
gex_matrix_schadendorf <- as.data.frame(gex_matrix_schadendorf)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_schadendorf, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "schadendorf_karlsson_9.csv", row.names=FALSE) # change name accordingly


# Gide ________________________________________________________________
gex_matrix_gide <- gex_matrix_gide[row.names(gex_matrix_gide) %in% common_genes_list,]
gex_matrix_gide <- as.data.frame(gex_matrix_gide)

prolif.karlsson <- get_signature_value_for_matrix(gex_matrix_gide, prolif_karlsson_in_common, max)
prolif.karlsson <- prolif.karlsson %>% dplyr::rename(prolif.karlsson = value)

write.csv(prolif.karlsson, "gide_karlsson_9.csv", row.names=FALSE) # change name accordingly


# After proliferation values from all tables have been combined
# divide each tumor type and data set into three proliferation groups 
combined$prolif.karlsson.group.intracancer <- "NA"
for (cancer.group in levels(combined$groups.to.analyze)) {
  print(cancer.group)
  # get only the data for that set and cancer type
  current.data <- subset(combined, groups.to.analyze == cancer.group)
  p_33 <- quantile(current.data$prolif.karlsson, .33)
  p_67 <- quantile(current.data$prolif.karlsson, .67)
  # if karl.value < 33th percentile, call it "Low", etc
  lows <- current.data[ which(current.data$prolif.karlsson < p_33), ] %>% pull(sample.id)
  intermediates <- current.data[ which(current.data$prolif.karlsson >= p_33 &
                                         current.data$prolif.karlsson <= p_67), ] %>% pull(sample.id)
  highs <- current.data[ which(current.data$prolif.karlsson > p_67), ] %>% pull(sample.id)
  # add to original table
  combined <- combined %>% mutate(prolif.karlsson.group.intracancer =
                                    case_when(sample.id %in% lows ~ "Low",
                                              sample.id %in% intermediates ~ "Intermediate",
                                              sample.id %in% highs ~ "High",
                                              TRUE ~ prolif.karlsson.group.intracancer))
}
combined$prolif.karlsson.group.intracancer <- factor(combined$prolif.karlsson.group.intracancer, 
                                                     levels=c("Low", "Intermediate", "High"))


# From now on, working with one master table containing every information
# e.g. CLAMS classification, proliferation value, overall survival time, molecular subtypes...


# Plot CLAMS proportion (Figure 1) -----
# base
combined %>% 
  group_by(groups.to.analyze, clams.class) %>% 
  tally() %>%
  pivot_wider(id_cols = groups.to.analyze, names_from = clams.class, values_from = n) %>% 
  mutate_at("TRU", ~replace(., is.na(.), 0))  %>% 
  mutate(tru_perc = (TRU / (TRU+NonTRU))) %>% mutate(total_samples = (TRU+NonTRU)) %>%
  mutate(tru_over_total = paste0('(', TRU, '/', total_samples, ')')) %>%
  arrange(desc(tru_perc)) %>% pivot_longer(cols=c("NonTRU", "TRU"), names_to = "clams.class") %>% 
  ggplot(aes(x=reorder(groups.to.analyze, tru_perc), y=value)) +
  geom_bar(position="fill", stat="identity", aes(fill=clams.class)) +
  theme_minimal() +
  coord_flip(clip = "off") +
  labs(x ="", y = "Samples") +
  geom_text(aes(x=groups.to.analyze, y=1.06, label=tru_over_total), size=3, color="grey50") +
  geom_text(aes(x=40, y=1.13, label="= n"), size=3, color="grey50") +
  geom_text(aes(x=41.2, y=0.05, label="TRU"), size=4, color="#0fbc6e") +
  geom_text(aes(x=41.2, y=0.9, label="NonTRU"), size=4, color="#0057b0") +
  theme(axis.text.x = element_text(size = 11),
        panel.grid.major.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        legend.position = "none") +
  scale_fill_manual(values=c("#0057b0", "#0fbc6e"), guide=guide_legend(reverse = TRUE), name="CLAMS") +
  scale_x_discrete(labels = full_names) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.25), labels=scales::percent)


# Overall survival analysis (Figure 2) -----
library(survival)
library(survminer)

combined <- convert_time(combined) # have time in both months and years
combined <- censor_data(combined, 5, "y") # censor data at 5 years
time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

os_sets <- c("THCA_TCGA", "KIP_TCGA", "PRAD_TCGA", "KICh_TCGA", "LUAD_TCGA",
             "KICC_TCGA", "PCPG_TCGA", "LGG_TCGA", "BRCA_SCAN-B", "ACC_TCGA",
             "BRCA_GOBO", "LIHC_TCGA", "PAAD_TCGA", "BRCA_TCGA")

for (os_set in os_sets) {
  print(os_set)
  type_subset <- subset(combined, groups.to.analyze == os_set)
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  # separating estimates by TRU or NonTRU
  fit_clams <- survfit(overall_surv_object~clams.class, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column]], 
                      type_subset[[event_column]])~type_subset$clams.class,))
  # simple plot for comparing info and making sure names are correct
  #print(ggsurvplot(fit_clams, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # better looking plot
  current_plot <- ggsurvplot(fit_clams, data=type_subset,
                             palette=c("#0fbc6e", "#0057b0"),
                             title=paste0("OS (", os_set, ")"),
                             xlab=paste0("Time (", time_type, ")"),
                             censor.shape=124, censor.size=3,
                             pval=TRUE, pval.coord=c(0,0.1),
                             #surv.median.line="hv",
                             risk.table=TRUE,
                             risk.table.fontsize = 4,
                             tables.theme = theme_survminer(font.main = 14),
                             legend="none", legend.title="CLAMS",
                             legend.labs=c("TRU", "NonTRU"))
  print(current_plot)
}

# correct p-values with Benjamini-Hochberg method
os_values <- read_csv('os_clams_values.csv')
os_values <- mutate(os_values, os.p.sign = ifelse(os_values$os.p < 0.05, "significant", "n"))
adjusted.bh <- round(p.adjust(os_values$os.p, method="fdr"), 4)
os_values <- cbind(os_values, adjusted.bh)
os_values <- mutate(os_values, adjusted.bh.sign = ifelse(os_values$adjusted.bh < 0.05, "significant", "n"))


# Cox proportional hazards analysis  -----
cox_sets <- c("KIP_TCGA", "LUAD_TCGA", "LGG_TCGA", "BRCA_SCAN-B", "BRCA_GOBO", "LIHC_TCGA", "BRCA_TCGA")

# Univariate
for (cox_set in cox_sets) {
  print(cox_set)
  type_subset <- subset(combined, groups.to.analyze == cox_set)
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~clams.class, data = type_subset)))
  forest <- ggforest(coxph(formula=overall_surv_object~clams.class, data=type_subset), 
                     main = paste("Hazard ratio", cox_set))
  print(forest)
}

# Multivariate (done for all groups in cox_sets)
multi_cox <- combined %>% filter(groups.to.analyze %in% cox_sets)

# Example: KIP (including gender, age, stage)
type_subset <- multi_cox %>% filter(groups.to.analyze == "KIP_TCGA")
overall_surv_object <- Surv(time=type_subset[[time_column]], 
                            event=type_subset[[event_column]])
print(summary(coxph(formula = overall_surv_object~clams.class+gender+Age+Stage.red, 
                    data = type_subset)))
forest <- ggforest(coxph(formula=overall_surv_object~clams.class+gender+Age+Stage.red, 
                         data=type_subset), 
                   main = "OS HR KIP_TCGA")
print(forest)


# Plot Proliferation boxplots (Figure 3) -----
# base
prolif_summary_tru <- combined[combined$clams.class == 'TRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_tru <- prolif_summary_tru %>% dplyr::rename(tru = n)
prolif_summary_nontru <- combined[combined$clams.class == 'NonTRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_nontru <- prolif_summary_nontru %>% dplyr::rename(nontru = n)
prolif_summary_tru_non <- merge(x=prolif_summary_nontru, y=prolif_summary_tru,
                                by="groups.to.analyze", all=TRUE)
rm(prolif_summary_nontru)
rm(prolif_summary_tru)

combined <- group_by(combined, groups.to.analyze, clams.class) %>% mutate(n=n())
above_10 <- subset(combined, n>=10)
tru_under_10 <- subset(combined, n<10)

max_tru <- combined %>%
  filter(clams.class == 'TRU') %>%
  select(prolif.karlsson) %>%
  max()
min_tru <- combined %>%
  filter(clams.class == 'TRU') %>%
  select(prolif.karlsson) %>%
  min()
max_tru_luad <- combined %>%
  filter(cancer.type == 'LUAD', clams.class == 'TRU') %>%
  select(prolif.karlsson) %>%
  max()
min_tru_luad <- combined %>%
  filter(cancer.type == 'LUAD', clams.class == 'TRU') %>%
  select(prolif.karlsson) %>%
  min()

combined %>%
  ggplot(aes(x=reorder(groups.to.analyze, prolif.karlsson), y=prolif.karlsson, fill=clams.class)) +
  geom_boxplot(data=above_10, position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#0fbc6e", "#0057b0"), name="CLAMS") +
  labs(data=combined, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_clams) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(breaks = c(2.6e+05, 4.8e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=6.2e+05, label=prolif_summary_tru_non$nontru, size=3, color="#0057b0", angle=90, hjust=0) +
  annotate("text", x=40.5, y=6.35e+05, label="= n", size=3, color="#0057b0", hjust=0) +
  annotate("text", x=0.75, y=6.7e+05, label="NonTRU", size=3, color="#0057b0", hjust=0) +
  annotate("text", x=0.75, y=1.3e+05, label="TRU", size=3, color="#0fbc6e", hjust=0) +
  annotate("text", x=36.5, y=1.55e+05, label="= n", size=3, color="#0fbc6e", hjust=0) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=1.5e+05, label=prolif_summary_tru_non$tru, size=3, color="#0fbc6e", angle=90, hjust=0) +
  geom_point(data=tru_under_10, size = 0.4, colour="#0fbc6e") +
  geom_hline(yintercept=max_luad_tru, linetype='dashed', color='gray', size=0.3) +
  geom_hline(yintercept=min_luad_tru, linetype='dashed', color='gray', size=0.3) +
  geom_hline(yintercept=max_tru, linetype='dotdash', color='#0fbc6e', size=0.3) +
  geom_hline(yintercept=min_tru, linetype='dotdash', color='#0fbc6e', size=0.3)



# Overall survival New proliferation  -----
prolif_os_groups <- combined %>% filter(treat_pred == 'no') %>% pull(groups.to.analyze) %>% unique()
prolif_os_groups <- sort(factor(prolif_os_groups))

list_plot <- vector(mode = "list", length = length(prolif_os_groups))
names(list_plot) <- paste0("plot_", levels(prolif_os_groups))
i <- 1

for (cancer_group in prolif_os_groups) {
  print(cancer_group)
  current.data <- combined %>% filter(groups.to.analyze == cancer_group)
  overall_surv_object <- Surv(time=current.data[[time_column]], 
                              event=current.data[[event_column]])
  # separating estimates by proliferation classes
  group_fit <- survfit(overall_surv_object~prolif.karlsson.group.intracancer, 
                       data=current.data)
  # testing
  print(survdiff(Surv(current.data[[time_column]], 
                      current.data[[event_column]])~current.data$prolif.karlsson.group.intracancer,))
  #print(ggsurvplot(tcga_fit, data=current.data, title=cancer_group), risk.table=TRUE))
  list_plot[[i]] <- ggsurvplot(group_fit, data=current.data,
                               palette=c("#fff533", "#ff875e", "#bc5090"),
                               title=cancer_group, 
                               xlab=paste0("Time (", time_type, ")"),
                               censor.shape=124, censor.size=3,
                               pval=TRUE, pval.coord=c(0,0.1),
                               #surv.median.line="hv",
                               risk.table=TRUE,
                               risk.table.fontsize = 4,
                               tables.theme = theme_survminer(font.main = 14),
                               legend="none", legend.title="Proliferation",
                               legend.labs=c("Low", "Intermediate", "High"))
  i <- i + 1
}

# base for Supplementary Material Figure 1
plots_1 <- arrange_ggsurvplots(list_plot[1:12], ncol = 3, nrow = 4, print=FALSE)
ggsave('prolif_OS_1.pdf', plots_1, width=20, height=28)
plots_2 <- arrange_ggsurvplots(list_plot[13:24], ncol = 3, nrow = 4, print=FALSE)
ggsave('prolif_OS_2.pdf', plots_2, width=20, height=28)
plots_3 <- arrange_ggsurvplots(list_plot[25:34], ncol = 3, nrow = 4, print=FALSE)
ggsave('prolif_OS_3.pdf', plots_3, width=20, height=28)


# Cox proportional hazards New proliferation  -----
no_treat <- combined %>% filter(treat_pred == 'no')
no_treat$prolif.karlsson.group.intracancer <- factor(no_treat$prolif.karlsson.group.intracancer,
                                                     levels = c('Low', 'Intermediate', 'High'))

for (type_to_analyze in levels(no_treat$groups.to.analyze)) {
  print(type_to_analyze)
  type_subset <- subset(no_treat, groups.to.analyze == type_to_analyze)
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~prolif.karlsson.group.intracancer, data = type_subset)))
  forest <- ggforest(coxph(formula=overall_surv_object~prolif.karlsson.group.intracancer, data=type_subset), 
                     main = paste("Hazard ratio", type_to_analyze))
  print(forest)
}

# correct p-values for High proliferation with Benjamini-Hochberg method
cox_5y_newprolif_table <- read_csv2("cox_5y_newprolif_table.csv")
high_prolif <- cox_5y_newprolif_table %>% filter(category == 'High') 
adjusted.bh <- round(p.adjust(high_prolif$p, method="fdr"), 4)
high_prolif <- cbind(high_prolif, adjusted.bh)
high_prolif <- mutate(high_prolif, adjusted.bh.sign = ifelse(high_prolif$adjusted.bh < 0.05, "significant", "n"))


# Plot New proliferation boxplots and forest plot (Figure 4)  -----
combined %>%
  ggplot(aes(x=reorder(groups.to.analyze, prolif.karlsson), y=prolif.karlsson, fill=prolif.karlsson.group.intracancer)) +
  geom_boxplot(data=combined, position="identity", outlier.shape=NA) +
  theme_minimal() +
  scale_fill_manual(values=c("#bc5090", "#ff875e","#fff533"), name="Proliferation") +
  labs(data=combined, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_cox_prolif_highlow) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 270, hjust=0), legend.position = "none",
        axis.text.y = element_text(angle = 270, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(breaks = c(2.6e+05, 4.8e+05), labels = c('Lower', 'Higher'))

library(forestplot)
cox_5y_newprolif_table <- read.csv("cox_5y_newprolif_table.csv")
cox_5y_newprolif_table <- cox_5y_newprolif_table %>% mutate(across(c(exp,lower_95,upper_95),log,.names="{col}.log"))
cox_5y_newprolif_table <- cox_5y_newprolif_table %>% mutate(group.names = paste(dataset, category))

to_plot <- cox_5y_newprolif_table[cox_5y_newprolif_table$significant == 'yes',]

forestplot(labeltext=to_plot$group.names, 
           mean=to_plot$exp.log,
           lower=to_plot$lower_95.log,
           upper=to_plot$upper_95.log,
           title='HR Proliferation',
           zero=0,
           boxsize=0.3,
           vertices = TRUE,
           xlab="HR and 95% CI",
           col=fpColors(box=c("darkblue", "darkred")),
           lwd.ci=2)



# Treatment prediction  -----

# Tabchy GSE20271 ________________________________________________________________
gse20271_clams_v_response <- combined %>%
                                filter(dataset == 'Tabchy') %>%
                                select(clams.class, Response) %>%
                                table()
gse20271_clams_v_response
fisher.test(gse20271_clams_v_response)


# Iwamoto GSE22093 ________________________________________________________________
gse22093_clams_v_response <- combined %>%
                              filter(dataset == 'Iwamoto') %>%
                              select(clams.class, Response) %>%
                              table()   
gse22093_clams_v_response
# can't do Fisher's


# Hess ________________________________________________________________
hess_clams_v_response <- combined %>%
                            filter(dataset == 'Hess') %>%
                            select(clams.class, Response) %>%
                            table()  
hess_clams_v_response
fisher.test(hess_clams_v_response)


# Riaz BMS038 ________________________________________________________________
riaz_clams_v_response <- combined %>%
                            filter(dataset == 'Riaz') %>%
                            select(clams.class, Response) %>%
                            table()  
# can't do Fisher's


# Liu ________________________________________________________________
schadendorf_clams_v_response <- combined %>%
                                    filter(dataset == 'Liu') %>%
                                    filter(Response != 'MR') %>%
                                    select(clams.class, Response) %>%
                                    table() 
schadendorf_clams_v_response
fisher.test(schadendorf_clams_v_response)


# Gide ________________________________________________________________
gide_clams_v_response <- combined %>%
                            filter(dataset == 'Gide') %>%
                            select(clams.class, Response) %>%
                            table()  
fisher.test(gide_clams_v_response)


# Supplementary Material Figure 2A -----
# base
subtypes_prolif <- combined %>% 
  filter(treat_pred == 'no' & cancer.type %in% c('BRCA', 'LGG', 'LIHC', 'KIP')) %>%
  pivot_longer(cols = c(PAM50, lgg.subtype, lihc.iCluster, lihc.hs, chen_taxonomy), 
               names_to = 'subtype.system', values_to = 'subtype', values_drop_na = TRUE) %>%
  select(dataset, cancer.type,  subtype.system, subtype, prolif.karlsson, groups.to.analyze) %>% 
  filter(subtype != 'Unclassified') %>%
  mutate(sub_group = paste(dataset, subtype, sep = '_'))
sample_size <- subtypes_prolif %>% group_by(sub_group) %>% summarize(num=n())
subtypes_prolif %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(subtype, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=prolif.karlsson)) +
  geom_violin() +
  theme_minimal() +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.y = element_text(angle = 90, hjust=0.5), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(2.5e+05, 5.5e+05), labels = c('Lower', 'Higher')) +
  xlab("Subtypes") +
  facet_wrap(~groups.to.analyze,nrow = 2,scales = "free_x")


# Supplementary Material Figure 2B -----
# BRCA Luminal A proliferation boxplots
sample_size <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'TCGA') %>%
  select(dataset, PAM50, prolif.karlsson, clams.class) %>% 
  group_by(clams.class) %>% 
  summarize(num=n())
p1 <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'TCGA') %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clams.class, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=prolif.karlsson, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#0fbc6e", "#0057b0"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  ggtitle("TCGA") +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10),
        plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = c(2.5e+05, 4.5e+05), labels = c('Lower', 'Higher')) +
  coord_cartesian(ylim = c(2.1e+05, 5e+05))

sample_size <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'GOBO') %>%
  select(dataset, PAM50, prolif.karlsson, clams.class) %>% 
  group_by(clams.class) %>% 
  summarize(num=n())
p2 <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'GOBO') %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clams.class, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=prolif.karlsson, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#0fbc6e", "#0057b0"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  ggtitle("GOBO") +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10),
        plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = c(2.5e+05, 4.5e+05), labels = c('Lower', 'Higher')) +
  coord_cartesian(ylim = c(2.1e+05, 5e+05))

sample_size <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'SCAN-B') %>%
  select(dataset, PAM50, prolif.karlsson, clams.class) %>% 
  group_by(clams.class) %>% 
  summarize(num=n())
p3 <- brca %>% 
  filter(PAM50 == 'LumA' & dataset == 'SCAN-B') %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clams.class, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=prolif.karlsson, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#0fbc6e", "#0057b0"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  ggtitle("SCAN-B") +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10),
        plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = c(2.5e+05, 4.5e+05), labels = c('Lower', 'Higher')) +
  coord_cartesian(ylim = c(2.1e+05, 5e+05))

sample_size <- brca %>% 
  filter(PAM50 == 'LumA') %>%
  select(dataset, PAM50, prolif.karlsson, clams.class) %>% 
  group_by(clams.class) %>% 
  summarize(num=n())
p4 <- brca %>% 
  filter(PAM50 == 'LumA') %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clams.class, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=prolif.karlsson, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#0fbc6e", "#0057b0"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  ggtitle("Three sets") +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=40,r=20,b=5,l=10),
        plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = c(2.5e+05, 4.5e+05), labels = c('Lower', 'Higher')) +
  coord_cartesian(ylim = c(2.1e+05, 5e+05))

plot_grid(p1, p2, p3, p4, ncol=4)