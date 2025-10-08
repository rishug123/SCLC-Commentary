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

# # Create cesa object that will include all data
cesa_all <- CESAnalysis(refset = "ces.refset.hg19")
# save_cesa(cesa_all, "cesa_all.rds")
# # cesa_all <- readRDS("cesa_all.rds")

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

# Summarize the amount of clinical and genomic data loaded in analysis

data_summary <- summarize_maf_sources(cesa_all)


# ---------------Assign COSMIC Mutational Signatures---------------------

# Define signature exclusions and use reference-based assignment (either "MutationalPatterns" or "deconstructSigs") to assign mutational signatures
assign_mutational_signatures("MutationalPatterns")

biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended==FALSE]
biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended==FALSE]

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

# --------Use dNdScv to find gene-level baseline mutation rates----------
# ---------------and estimate gene-level cancer effect-------------------


#Calculate gene_mutation_rates
cesa_paper <- gene_mutation_rates(cesa_paper, covariates = ces.refset.hg19$covariates$lung)

cesa_external <- gene_mutation_rates(cesa_external, covariates = ces.refset.hg19$covariates$lung)

# Store top dNdScv genes, sorted by significance
dndscv_results_paper <- cesa_paper$dNdScv_results$rate_grp_1
sig_genes_paper <- dndscv_results_paper[dndscv_results_paper$qallsubs_cv < 0.05, ]

dndscv_results_external <- cesa_external$dNdScv_results$rate_grp_1
sig_genes_external <- dndscv_results_external[dndscv_results_external$qallsubs_cv < 0.05, ]

#--------Find the variant cancer-effect of defined variant lists--------

# The following will include every single variant with a MAF frequency >1
cesa_external <- ces_variant(cesa = cesa_external, run_name = "recurrent_external")
cesa_paper <- ces_variant(cesa = cesa_paper, run_name = "recurrent_paper")

## -------------------- Effect sizes attributed to signatures -----

signature_attribution_external <- mutational_signature_effects(cesa = cesa_external, effects = cesa_external$selection$recurrent_external)

signature_attribution_paper <- mutational_signature_effects(cesa = cesa_paper, effects = cesa_paper$selection$recurrent_paper)

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

# Tried just using Li et al., 2025 signature data but could not replicate the survival findings

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
