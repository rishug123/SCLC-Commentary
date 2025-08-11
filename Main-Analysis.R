setwd("..") # the location the r script will save

# -----Define reference data set and load data into cesa----------


# Create cesa object that will include external data 
cesa_external <- CESAnalysis(refset = "ces.refset.hg19")
save_cesa(cesa_external, "cesa_external.rds")
# cesa_external <- readRDS("cesa_external.rds")

# Create cesa object that will include only data in the paper
cesa_paper <- CESAnalysis(refset = "ces.refset.hg19")
save_cesa(cesa_paper, "cesa_paper_new.rds")
# cesa_paper <- readRDS("cesa_paper.rds")

# Create cesa object that will include all data
cesa_all <- CESAnalysis(refset = "ces.refset.hg19")
save_cesa(cesa_all, "cesa_all.rds")
# cesa_all <- readRDS("cesa_all.rds")

# Check for duplicate data
maf_list <- list(dt1 = george_et_al_maf, dt2 = rudin_et_al_maf, dt3 = jiang_et_al_maf, dt4 = zhou_et_al_maf)
maf_overlap <- check_sample_overlap(maf_list)
nrow(maf_overlap[variants_shared > 2])

# Load WXS data
cesa_all <- load_maf(cesa = cesa_all, maf = george_et_al_maf, maf_name = "george_et_al")
cesa_paper <- load_maf(cesa = cesa_paper, maf = george_et_al_maf, maf_name = "george_et_al")

cesa_all <- load_sample_data(cesa_all, george_et_al_clinical)
cesa_paper <- load_sample_data(cesa_paper, george_et_al_clinical)

cesa_all <- load_maf(cesa = cesa_all, maf = rudin_et_al_maf, maf_name = "rudin_et_al")
cesa_paper <- load_maf(cesa = cesa_paper, maf = rudin_et_al_maf, maf_name = "rudin_et_al")

# External data
cesa_external <- load_maf(cesa = cesa_external, maf = zhou_et_al_maf, maf_name = "zhou_et_al")
cesa_external <- load_sample_data(cesa_external, zhou_et_al_clinical)

cesa_external <- load_maf(cesa = cesa_external, maf = jiang_et_al_maf, maf_name = "jiang_et_al")
cesa_external <- load_sample_data(cesa_external, jiang_et_al_clinical)

cesa_all <- load_maf(cesa = cesa_all, maf = zhou_et_al_maf, maf_name = "zhou_et_al")
cesa_all <- load_sample_data(cesa_all, zhou_et_al_clinical)

cesa_all <- load_maf(cesa = cesa_all, maf = jiang_et_al_maf, maf_name = "jiang_et_al")
cesa_all <- load_sample_data(cesa_all, jiang_et_al_clinical)

# # load paper targeted data set into the cesa object
# 
#   # paper
# 
# cesa_paper <- load_maf(cesa_paper, maf = IMPACT468_msk_met_maf, coverage = "targeted", covered_regions = IMPACT468_msk_met_coverage, covered_regions_name = "IMPACT468", maf_name = "msk_met")
# cesa_paper <- load_maf(cesa_paper, maf = IMPACT410_msk_met_maf, coverage = "targeted", covered_regions = IMPACT468_msk_met_coverage, covered_regions_name = "IMPACT410", maf_name = "msk_met")
# cesa_paper <- load_maf(cesa_paper, maf = IMPACT341_msk_met_maf, coverage = "targeted", covered_regions = IMPACT341_msk_met_coverage, covered_regions_name = "IMPACT341", maf_name = "msk_met")
# 
# cesa_paper <- load_sample_data(cesa_paper, msk_met_clinical_txt_file)
# 
#   # all
# 
# cesa_all <- load_maf(cesa_all, maf = IMPACT468_msk_met_maf, coverage = "targeted", covered_regions = IMPACT468_msk_met_coverage, covered_regions_name = "IMPACT468", maf_name = "msk_met")
# cesa_all <- load_maf(cesa_all, maf = IMPACT410_msk_met_maf, coverage = "targeted", covered_regions = IMPACT468_msk_met_coverage, covered_regions_name = "IMPACT410", maf_name = "msk_met")
# cesa_all <- load_maf(cesa_all, maf = IMPACT341_msk_met_maf, coverage = "targeted", covered_regions = IMPACT341_msk_met_coverage, covered_regions_name = "IMPACT341", maf_name = "msk_met")
# 
# cesa_all <- load_sample_data(cesa_all, msk_met_clinical_txt_file)

# Define the samples used in Li et al., 2025 (Nature)
genes_identified_in_paper <- c( "KMT2C", "THSD7B", "WDR87", "UNC13A", "SPATA31D1", "SLC4A10", "PKD1L1", "ANKS1B", "LTN1", "RSF1", "CCT8L2", "VPS13D", "USP13", "NAALAD2", "NDST3", "ACACB", "POM121L12", "CEP128", "BNC2", "DNAH10", "NTRK2", "DOK6", "SAGE1", "OR2T27", "ST6GALNAC3", "FMR1", "TBX18", "UGGT2", "FAM47B", "SPRY3", "CECR2", "ADAMTS18", "OR14I1", "GPR126", "SLC5A11", "SEMA5A", "TPH2", "FOLH1", "TP53")

# Summarize the amount of clinical and genomic data loaded in analysis

data_summary <- summarize_maf_sources(cesa_all)


# ---------------Assign COSMIC Mutational Signatures---------------------

# Define signature exclusions and use reference-based assignment (either "MutationalPatterns" or "deconstructSigs") to assign mutational signatures
assign_mutational_signatures("MutationalPatterns")

biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended==FALSE]
biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended==FALSE]
biological_weights_all <- cesa_all$mutational_signatures$biological_weights[group_avg_blended==FALSE]

# non-bootstrapped optimal cutoff for paper data

# Load clinical data and prepare for analysis
clinical_data_paper <- cesa_paper$samples %>%
  filter(
    !is.na(`Status.(at.time.of.last.follow-up)`),
    !is.na(`overall_survival.(months)`)
  ) %>%
  mutate(
    status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)
  ) %>%
  select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)

# Merge mutation and clinical data
paper_optimal_cutoff <- biological_weights_paper %>%
  # select(total_snvs, Unique_Patient_Identifier, SBS4, SBS13, SBS5, SBS29, SBS3,SBS24, SBS87) %>%
  inner_join(clinical_data_paper, by = "Unique_Patient_Identifier")

# Perform optimal cutoff tests (non-bootstrap)
SBS4_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS4,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)

SBS4_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS4,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)

SBS13_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS13,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)

SBS13_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS13,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)

SBS3_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS3,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)
SBS3_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS3, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
)
SBS87_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS87,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)
SBS87_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS87, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
)

SBS24_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS24,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)
SBS24_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS24, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
)

SBS31_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(`overall_survival.(months)`, status_numeric) ~ SBS31,
  data = paper_optimal_cutoff,
  smethod = "LogRank", pmethod = "HL"
)
SBS31_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
  Surv(total_snvs) ~ SBS31, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
)

# --------Use dNdScv to find gene-level baseline mutation rates----------
# ---------------and estimate gene-level cancer effect-------------------


#Calculate gene_mutation_rates
cesa_paper <- gene_mutation_rates(cesa_paper, covariates = ces.refset.hg19$covariates$lung)

cesa_external <- gene_mutation_rates(cesa_external, covariates = ces.refset.hg19$covariates$lung)

cesa_all <- gene_mutation_rates(cesa_all, covariates = ces.refset.hg19$covariates$lung)

# Store top dNdScv genes, sorted by significance
dndscv_results_paper <- cesa_paper$dNdScv_results$rate_grp_1
sig_genes_paper <- dndscv_results_paper[dndscv_results_paper$qallsubs_cv < 0.05, ]

dndscv_results_external <- cesa_external$dNdScv_results$rate_grp_1
sig_genes_external <- dndscv_results_external[dndscv_results_external$qallsubs_cv < 0.05, ]

dndscv_results_all <- cesa_all$dNdScv_results$rate_grp_1
sig_genes_all <- dndscv_results_all[dndscv_results_all$qallsubs_cv < 0.05, ]



# #---------Defining variants to calculate cancer effect--------------------
# 
# # Significant genes, identified by gene-level effect size calculated by dndscv
# dndscv_estimated_variants_external <- cesa_external$variants[cesa_external$variants$gene %in% sig_genes_external$gene_name,]
# dndscv_estimated_variants_paper <- cesa_paper$variants[cesa_paper$variants$gene %in% sig_genes_paper$gene_name,]
# dndscv_estimated_variants_all <- cesa_all$variants[cesa_all$variants$gene %in% sig_genes_all$gene_name,]
# 
# # Identified prognostic genes in the paper
# tested_variants_external <- cesa_external$variants[cesa_external$variants$gene %in% genes_identified_in_paper,]
# tested_variants_paper <- cesa_paper$variants[cesa_paper$variants$gene %in% genes_identified_in_paper,]
# tested_variants_all <- cesa_all$variants[cesa_all$variants$gene %in% genes_identified_in_paper,]
# 

# # ---------Use trinucleotide and gene level mutation estimates----------- 
# #-------------------as well as lung covariate information ----------------
# #------------------to calculate baseline mutation rate-------------------
# 
# # A few random samples in general
# samples_to_check_external <- sample(cesa_external$samples$Unique_Patient_Identifier, 3)
# samples_to_check_paper <- sample(cesa_paper$samples$Unique_Patient_Identifier, 3)
# samples_to_check_all <- sample(cesa_all$samples$Unique_Patient_Identifier, 3)
# 
# baseline_mutation_rates(cesa = cesa_external, variant_ids = dndscv_estimated_variants_external$variant_id, samples = samples_to_check_external)
# baseline_mutation_rates(cesa = cesa_paper, variant_ids = dndscv_estimated_variants_paper$variant_id, samples = samples_to_check_paper)
# baseline_mutation_rates(cesa = cesa_all, variant_ids = dndscv_estimated_variants_all$variant_id, samples = samples_to_check_all)
# 
# baseline_mutation_rates(cesa = cesa_external, variant_ids = tested_variants_external$variant_id, samples = samples_to_check_external)
# baseline_mutation_rates(cesa = cesa_paper, variant_ids = tested_variants_paper$variant_id, samples = samples_to_check_paper)
# baseline_mutation_rates(cesa = cesa_all, variant_ids = tested_variants_all$variant_id, samples = samples_to_check_all)


#--------Find the variant cancer-effect of defined variant lists--------

# The following will include every single variant with a MAF frequency >1
cesa_external <- ces_variant(cesa = cesa_external, run_name = "recurrent_external")
cesa_paper <- ces_variant(cesa = cesa_paper, run_name = "recurrent_paper")
cesa_all <- ces_variant(cesa = cesa_all, run_name = "recurrent_all")

# # Testing to see if any of the top 100 effect size genes from recurrent analysis have prognostic significance in cesa_external
# recurrent_external_genes_data <- get_prognostic_significance(cesa_external$selection$recurrent_external, top_x_genes, cox_regression_df, filter_by = "P")
# recurrent_external_genes_data$summary
# 
# # Testing to see if any of the top 100 effect size genes from recurrent analysis have prognostic significance in cesa_all
# recurrent_all_genes_data <- get_prognostic_significance(cesa_all$selection$recurrent_all, top_x_genes, cox_regression_df, filter_by = "P")
# recurrent_all_genes_data$summary
# 
# # Testing to see if any of the top 100 genes from recurrent analysis have prognostic significance in cesa_paper
# recurrent_paper_genes_data <- get_prognostic_significance(cesa_paper$selection$recurrent_paper, top_x_genes, cox_regression_df, filter_by = "P")
# recurrent_paper_genes_data$summary
# 
# # Variant cancer-effect for dndscv identified genes
# 
# cesa_external <- ces_variant(cesa = cesa_external, run_name = "dndscv_external", variants = dndscv_estimated_variants_external)
# cesa_paper <- ces_variant(cesa = cesa_paper, run_name = "dndscv_paper", variants = dndscv_estimated_variants_paper)
# cesa_all <- ces_variant(cesa = cesa_all, run_name = "dndscv_all", variants = dndscv_estimated_variants_all)
# 
# # Testing to see if any of the top 1701 effect size genes from dndscv analysis have prognostic significance in cesa_external
# dndscv_external_genes_data <- get_prognostic_significance(cesa_external$selection$dndscv_external, nrow(sig_genes_external), cox_regression_df, filter_by = "P")
# dndscv_external_genes_data$summary
# # dndscv_external_genes_data$result %>%
# #   arrange((multivariate_fdr))
# 
# # Testing to see if any of the top effect size genes from dndscv analysis have prognostic significance in cesa_all
# dndscv_all_genes_data <- get_prognostic_significance(cesa_all$selection$dndscv_all, nrow(sig_genes_all), cox_regression_df, filter_by = "P")
# dndscv_all_genes_data$summary
# # dndscv_all_genes_data$result %>%
# #   arrange((multivariate_fdr))
# 
# # Testing to see if any of the top 4 effect size genes from dndscv analysis have prognostic significance in cesa_paper
# dndscv_paper_genes_data <- get_prognostic_significance(cesa_paper$selection$dndscv_paper, nrow(sig_genes_paper), cox_regression_df, filter_by = "P")
# dndscv_paper_genes_data$summary
# 
# # Variant cancer effect for Li et al., 2025 identified prognostic genes
# cesa_external <- ces_variant(cesa = cesa_external, run_name = "tested_external", variants = tested_variants_external)
# cesa_paper <- ces_variant(cesa = cesa_paper, run_name = "tested_paper", variants = tested_variants_paper)
# cesa_all <- ces_variant(cesa = cesa_all, run_name = "tested_all", variants = tested_variants_all)
# 
# # Testing to see if any of the top 100 effect size genes from Li et al., 2025 identified prognostic gene analysis have prognostic significance in cesa_external
# tested_external_genes_data <- get_prognostic_significance(cesa_external$selection$tested_external, top_x_genes, cox_regression_df, filter_by = "P")
# tested_external_genes_data$summary
# 
# # Testing to see if any of the top 100 effect size genes from Li et al., 2025 identified prognostic gene analysis have prognostic significance in cesa_all
# tested_all_genes_data <- get_prognostic_significance(cesa_all$selection$tested_all, top_x_genes, cox_regression_df, filter_by = "P")
# tested_all_genes_data$summary
# 
# # Testing to see if any of the top 100 effect size genes from Li et al., 2025 identified prognostic gene analysis have prognostic significance in cesa_paper
# tested_paper_genes_data <- get_prognostic_significance(cesa_paper$selection$tested_paper, top_x_genes, cox_regression_df, filter_by = "P")
# tested_paper_genes_data$summary


## -------------------- Effect sizes attributed to signatures -----

signature_attribution_external <- mutational_signature_effects(cesa = cesa_external, effects = cesa_external$selection$recurrent_external)

avg_by_variant_SBS13_external <- get_gene_prognosis_from_signature_attribution(cesa_external, signature_attribution_external, "SBS13")
avg_by_variant_SBS4_external <- get_gene_prognosis_from_signature_attribution(cesa_external, signature_attribution_external, "SBS4")
avg_by_variant_SBS5_external <- get_gene_prognosis_from_signature_attribution(cesa_external, signature_attribution_external, "SBS5")


signature_attribution_paper <- mutational_signature_effects(cesa = cesa_paper, effects = cesa_paper$selection$recurrent_paper)

avg_by_variant_SBS13_paper <- get_gene_prognosis_from_signature_attribution(cesa_paper, signature_attribution_paper, "SBS13")
avg_by_variant_SBS4_paper <- get_gene_prognosis_from_signature_attribution(cesa_paper, signature_attribution_paper, "SBS4")
avg_by_variant_SBS5_paper <- get_gene_prognosis_from_signature_attribution(cesa_paper, signature_attribution_paper, "SBS5")

view_signature_attribution <- cesa_paper$selection$recurrent_paper %>%
  inner_join(signature_attribution_paper$mutational_sources$average_by_variant, by = "variant_id") %>%
  # select(variant_name, selection_intensity, SBS4, SBS13, gene) %>%
  filter(gene %in% c("PIK3CA","PITX2", "GALNT5", "PRSS57", "WDR65", "ADHGEF1", "NAALAD2", "NDST3", "GPR126", "FAM47B", "MGAT4C", "NTRK2"))

signature_attribution_paper$mutational_sources$average_by_variant$variant_id


signature_attribution_all <- mutational_signature_effects(cesa = cesa_all, effects = cesa_all$selection$recurrent_all)

avg_by_variant_SBS13_all <- get_gene_prognosis_from_signature_attribution(cesa_all, signature_attribution_all, "SBS13")
avg_by_variant_SBS4_all <- get_gene_prognosis_from_signature_attribution(cesa_all, signature_attribution_all, "SBS4")
avg_by_variant_SBS5_all <- get_gene_prognosis_from_signature_attribution(cesa_all, signature_attribution_all, "SBS5")

view_signature_attribution <- cesa_all$selection$recurrent_all %>%
  inner_join(signature_attribution_all$mutational_sources$average_by_variant, by = "variant_id") %>%
  select(variant_name, selection_intensity, SBS4, SBS13, SBS5, gene) %>%
  filter(gene %in% c("PIK3CA","PITX2", "GALNT5", "PRSS57", "WDR65", "ADHGEF1", "NAALAD2", "NDST3", "GPR126", "FAM47B", "MGAT4C", "NTRK2"))

sex_count <- cesa_paper$samples %>% filter(!is.na(`overall_survival.(months)`)) %>% filter(!is.na(age)) %>% filter(coverage == "exome")
sum(sex_count$sex == "female")
nrow(sex_count)


# #-------Estimate gene cancer-effect using compound variant analysis-------
# 
# # external
# 
# # Note: see Functions
# 
# # Set up compound variant analysis
# 
# external_cov = c(cesa_external$coverage_ranges$exome, cesa_external$coverage_ranges$targeted)
# external_cov = Reduce(GenomicRanges::intersect, external_cov)
# 
# 
# # Genes within the 13 ranges (only required if including GENIE database)
# # top_genes_for_compound_variant <- c("CTNNB1", "PIK3CA", "BRAF", "PTEN", "KRAS", "TP53")
# 
# # Use compound variant set to find out the gene level selection of top 35, recurrent genes
# 
# top_genes_for_compound_variant <- cesa_external$selection$recurrent_external %>%
#   arrange(desc(selection_intensity)) %>%
#   distinct(gene, .keep_all = FALSE) %>%
#   slice_head(n = 35) %>%
#   pull(gene)
# 
# # Define data.table with all variants associated with specified genes
# variants <- select_variants(cesa_external, genes = top_genes_for_compound_variant, gr = external_cov)
# 
# # Process all genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(top_genes_for_compound_variant, genes_info)
# 
# # # Apply the prevalence filter (>25) and compound variant definition as in your original code
# # comp_genes_high_prevalence <- for_comp %>%
# #   group_by(gene) %>%
# #   summarize(occurrence = sum(maf_prevalence)) %>%
# #   filter(occurrence >= 5)
# # 
# # for_comp <- for_comp %>%
# #   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(cesa = cesa_external, variant_table = for_comp, by = "gene", merge_distance = Inf)
# 
# cesa_external <- ces_variant(cesa = cesa_external, run_name = "gene_level_external", variants = compound)
# gene_level_external <- cesa_external$selection$gene_level_external
# colnames(gene_level_external)[colnames(gene_level_external) == "variant_name"] <- "gene"
# gene_level_external_data <- get_prognostic_significance(gene_level_external, nrow(gene_level_external), cox_regression_df, filter_by = "P")
# 
# gene_level_external_data$summary
# 
# 
# # Use compound variant analysis for paper identified genes
# 
# # Select variants
# 
# variants <- select_variants(cesa_external, genes = genes_identified_in_paper, gr = external_cov)
# 
# # Process external genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(
#   top_genes_for_compound_variant = genes_identified_in_paper,
#   genes_info = genes_info
# )
# 
# # Optional prevalence filter (uncomment if needed)
# comp_genes_high_prevalence <- for_comp %>%
#   group_by(gene) %>%
#   summarize(occurrence = sum(maf_prevalence)) %>%
#   filter(occurrence > 5)
# #
# for_comp <- for_comp %>%
#   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(
#   cesa = cesa_external,
#   variant_table = for_comp,
#   by = "gene",
#   merge_distance = Inf
# )
# 
# # Add compound variant annotations to the CESAnalysis object
# cesa_external <- ces_variant(
#   cesa = cesa_external,
#   run_name = "gene_level_external_for_paper",
#   variants = compound
# )
# 
# # all
# 
# # Set up compound variant analysis
# 
# all_cov = c(cesa_all$coverage_ranges$exome, cesa_all$coverage_ranges$targeted)
# all_cov = Reduce(GenomicRanges::intersect, all_cov)
# 
# # Use compound variant set to find out the gene level selection of top 35, recurrent genes
# 
# top_genes_for_compound_variant <- cesa_all$selection$recurrent_all %>%
#   arrange(desc(selection_intensity)) %>%
#   distinct(gene, .keep_all = FALSE) %>%
#   slice_head(n = 35) %>%
#   pull(gene)
# 
# # Define data.table with all variants associated with specified genes
# variants <- select_variants(cesa_all, genes = top_genes_for_compound_variant, gr = all_cov)
# 
# # Process all genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(top_genes_for_compound_variant, genes_info)
# 
# # # Apply the prevalence filter (>25) and compound variant definition as in your original code
# # comp_genes_high_prevalence <- for_comp %>%
# #   group_by(gene) %>%
# #   summarize(occurrence = sum(maf_prevalence)) %>%
# #   filter(occurrence >= 5)
# # 
# # for_comp <- for_comp %>%
# #   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(cesa = cesa_all, variant_table = for_comp, by = "gene", merge_distance = Inf)
# 
# cesa_all <- ces_variant(cesa = cesa_all, run_name = "gene_level_all", variants = compound)
# gene_level_all <- cesa_all$selection$gene_level_all
# colnames(gene_level_all)[colnames(gene_level_all) == "variant_name"] <- "gene"
# gene_level_all_data <- get_prognostic_significance(gene_level_all, nrow(gene_level_all), cox_regression_df, filter_by = "P")
# 
# gene_level_all_data$summary
# 
# 
# # Use compound variant analysis for paper identified genes
# 
# # Select variants
# 
# variants <- select_variants(cesa_all, genes = genes_identified_in_paper, gr = all_cov)
# 
# # Process external genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(
#   top_genes_for_compound_variant = genes_identified_in_paper,
#   genes_info = genes_info
# )
# 
# # Optional prevalence filter (uncomment if needed)
# comp_genes_high_prevalence <- for_comp %>%
#   group_by(gene) %>%
#   summarize(occurrence = sum(maf_prevalence)) %>%
#   filter(occurrence > 5)
# #
# for_comp <- for_comp %>%
#   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(
#   cesa = cesa_all,
#   variant_table = for_comp,
#   by = "gene",
#   merge_distance = Inf
# )
# 
# # Add compound variant annotations to the CESAnalysis object
# cesa_all <- ces_variant(
#   cesa = cesa_all,
#   run_name = "gene_level_all_for_paper",
#   variants = compound
# )
# 
# # gene_level_all_for_paper <- cesa_paper$selection$gene_level_all_for_paper
# # colnames(gene_level_all_for_paper)[colnames(gene_level_all_for_paper) == "variant_name"] <- "gene"
# # gene_level_all_for_paper_data <- get_prognostic_significance(gene_level_all_for_paper, nrow(gene_level_all_for_paper), cox_regression_df, filter_by = "P")
# # gene_level_all_for_paper_data$summary
# 
# # Paper
# 
# # Note: see Functions
# 
# # Set up compound variant analysis
# 
# paper_cov = c(cesa_paper$coverage_ranges$exome, cesa_paper$coverage_ranges$targeted)
# paper_cov = Reduce(GenomicRanges::intersect, paper_cov)
# 
# 
# # Use compound variant set to find out the gene level selection of top, recurrent genes
# 
# top_genes_for_compound_variant <- cesa_paper$selection$recurrent_paper %>%
#   filter(!is.na(gene)) %>%                            # Remove NA gene names
#   arrange(desc(selection_intensity)) %>%              # Sort by selection intensity
#   distinct(gene, .keep_all = FALSE) %>%               # Keep only unique genes
#   slice_head(n = 35) %>%                              # Select top 50
#   pull(gene)                                          # Extract as vector
# 
# 
# # Define data.table with all variants associated with specified genes
# variants <- select_variants(cesa_paper, genes = top_genes_for_compound_variant, gr = paper_cov)
# 
# # Process external genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(top_genes_for_compound_variant, genes_info)
# 
# # # Apply the prevalence filter (>25) and compound variant definition as in your original code
# comp_genes_high_prevalence <- for_comp %>%
#   group_by(gene) %>%
#   summarize(occurrence = sum(maf_prevalence)) %>%
#   filter(occurrence >= 5)
# 
# for_comp <- for_comp %>%
#   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(cesa = cesa_paper, variant_table = for_comp, by = "gene", merge_distance = Inf)
# 
# cesa_paper <- ces_variant(cesa = cesa_paper, run_name = "gene_level_paper", variants = compound)
# gene_level_paper_recurrent <- cesa_paper$selection$gene_level_paper
# colnames(gene_level_paper_recurrent)[colnames(gene_level_paper_recurrent) == "variant_name"] <- "gene"
# gene_level_paper_recurrent_data <- get_prognostic_significance(gene_level_paper_recurrent, nrow(gene_level_paper_recurrent), cox_regression_df, filter_by = "P")
# gene_level_paper_recurrent_data$summary
# 
# 
# # Use compound variant analysis for paper identified genes
# 
# 
# # Select variants
# variants <- select_variants(cesa_paper, genes = genes_identified_in_paper, gr = paper_cov)
# 
# # Process all genes to define a variant list for compound variant analysis
# for_comp <- process_all_genes(
#   top_genes_for_compound_variant = genes_identified_in_paper,
#   genes_info = genes_info
# )
# 
# # Optional prevalence filter (uncomment if needed)
# # comp_genes_high_prevalence <- for_comp %>%
# #   group_by(gene) %>%
# #   summarize(occurrence = sum(maf_prevalence)) %>%
# #   filter(occurrence > 5)
# # #
# # for_comp <- for_comp %>%
# #   filter(gene %in% comp_genes_high_prevalence$gene)
# 
# # Define compound variants
# compound <- define_compound_variants(
#   cesa = cesa_paper,
#   variant_table = for_comp,
#   by = "gene",
#   merge_distance = Inf
# )
# 
# # Add compound variant annotations to the CESAnalysis object
# cesa_paper <- ces_variant(
#   cesa = cesa_paper,
#   run_name = "gene_level_paper_for_paper",
#   variants = compound
# )
# colnames(cesa_paper$selection$gene_level_paper)
# gene_level_paper <- cesa_paper$selection$gene_level_paper
# colnames(gene_level_paper)[colnames(gene_level_paper) == "variant_name"] <- "gene"
# gene_level_paper_data <- get_prognostic_significance(gene_level_paper, nrow(gene_level_paper), cox_regression_df, filter_by = "P")
# gene_level_paper_data$summary
# 
# gene_level_paper_for_paper <- cesa_paper$selection$gene_level_paper_for_paper
# colnames(gene_level_paper_for_paper)[colnames(gene_level_paper_for_paper) == "variant_name"] <- "gene"
# gene_level_paper_for_paper_data <- get_prognostic_significance(gene_level_paper_for_paper, nrow(gene_level_paper_for_paper), cox_regression_df, filter_by = "P")
# gene_level_paper_for_paper_data$summary
# 
# 

#--------------Estimate prognosis linked to mutational signatures ------------

survival_data_external <- cesa_external$samples %>%
  select(
    Unique_Patient_Identifier, staging,
    overall_survival_months = `overall_survival.(months)`,
    chemo = `chemotherapy.(yes/no)`,
    sclc_treatment = previous.therapeutic.treatment.for.SCLC,
    status = `Status.(at.time.of.last.follow-up)`
  ) %>%
  mutate(
    status_numeric = ifelse(status == "dead", 1, 0)
  ) %>%
  filter(!is.na(overall_survival_months), !is.na(status_numeric))

survival_data_paper <- cesa_paper$samples %>%
  select(
    Unique_Patient_Identifier,staging,
    overall_survival_months = `overall_survival.(months)`,
    chemo = `chemotherapy.(yes/no)`,
    sclc_treatment = previous.therapeutic.treatment.for.SCLC,
    status = `Status.(at.time.of.last.follow-up)`
  ) %>%
  mutate(
    status_numeric = ifelse(status == "dead", 1, 0)
  ) %>%
  filter(!is.na(overall_survival_months), !is.na(status_numeric))

survival_data_all <- cesa_all$samples %>%
  select(
    Unique_Patient_Identifier,staging,
    overall_survival_months = `overall_survival.(months)`,
    chemo = `chemotherapy.(yes/no)`,
    sclc_treatment = previous.therapeutic.treatment.for.SCLC,
    status = `Status.(at.time.of.last.follow-up)`
  ) %>%
  mutate(
    status_numeric = ifelse(status == "dead", 1, 0)
  ) %>%
  filter(!is.na(overall_survival_months), !is.na(status_numeric))

# Complete the following analyses: Kaplan-Meier analysis using 0.25 binary, Kaplan-Meier analysis using optimal binary
SBS5_results <- analyze_signature_complete("SBS5")
SBS4_results <- analyze_signature_complete("SBS4")
SBS13_results <- analyze_signature_complete("SBS13")

# Hazard ratios w/ 95% CI
summary(SBS4_results$paper$binary$hazard_0.25)
summary(SBS4_results$paper$binary$hazard_optimal)
summary(SBS4_results$external$binary$hazard_0.25)
summary(SBS4_results$external$binary$hazard_optimal)
summary(SBS13_results$paper$binary$hazard_0.25)
summary(SBS13_results$paper$binary$hazard_optimal)
summary(SBS13_results$external$binary$hazard_0.25)
summary(SBS13_results$external$binary$hazard_optimal)

# tested association between signatures with high average proportion of TMB and overall survival

# Significant
SBS24_results <- analyze_signature_complete("SBS24", cutoff = round(SBS24_paper_optimal_TMB_cutoff_for_paper_non_bs$estimate, digits = 2))
summary(SBS24_results$paper$binary$hazard_optimal)

# Not significant
SBS3_results <- analyze_signature_complete("SBS3", cutoff = round(SBS3_paper_optimal_survival_cutoff_for_paper_non_bs$estimate, digits =2))
SBS87_results <- analyze_signature_complete("SBS87", cutoff = round(SBS87_paper_optimal_survival_cutoff_for_paper_non_bs$estimate, digits = 2))
SBS31_results <- analyze_signature_complete("SBS31", cutoff = round(SBS31_paper_optimal_survival_cutoff_for_paper_non_bs$estimate, digits = 2))

# Tried just using their signature data but could not replicate the survival findings

# Subset of paper_signatures in signature_analysis_paper
survival_data_paper$Unique_Patient_Identifier

paper_signatures <- paper_signatures[paper_signatures$Unique_Patient_Identifier %in% survival_data_paper$Unique_Patient_Identifier,]
paper_sig_survival_data_paper <- survival_data_paper[survival_data_paper$Unique_Patient_Identifier %in% paper_signatures$Unique_Patient_Identifier, ]

# Merge clinical and signature data
paper_signatures_optimal_cutoff <- paper_signatures %>%
  inner_join(paper_sig_survival_data_paper, by = "Unique_Patient_Identifier") %>%
  select(Unique_Patient_Identifier, SBS4, SBS5, SBS13, overall_survival_months, status_numeric)

# Compute optimal cutoff using maxstat.test
SBS4_paper_signatures_optimal_cutoff <- maxstat.test(
  Surv(overall_survival_months, status_numeric) ~ SBS4, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
  pmethod = "HL"
)
SBS13_paper_signatures_optimal_cutoff <- maxstat.test(
  Surv(overall_survival_months, status_numeric) ~ SBS13, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
  pmethod = "HL"
)

SBS5_paper_signatures_optimal_cutoff <- maxstat.test(
  Surv(overall_survival_months, status_numeric) ~ SBS5, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
  pmethod = "HL"
)

sig_result_SBS4 <- analyze_signature(signature_name = "SBS4", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)
summary(sig_result_SBS4$hazard_0.25)

sig_result_SBS13 <- analyze_signature(signature_name = "SBS13", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)
summary(sig_result_SBS13$hazard_0.25)

sig_result_SBS5 <- analyze_signature(signature_name = "SBS5", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)


#-----------------------Estimating prognosis for high effect genes--------------------
# # Placeholder for the get_prognostic_significance results
# 
# dndscv_external_genes_data
# dndscv_all_genes_data
# dndscv_paper_genes_data
# tested_external_genes_data 
# tested_all_genes_data
# tested_paper_genes_data 
# gene_level_external_data
# gene_level_all_data
# gene_level_paper_recurrent_data
# gene_level_paper_for_paper_data

