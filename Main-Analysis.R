setwd("..") # the location the r script will save

# -----Define reference data set and load data into cesa----------


# Create cesa object that will include external data 
cesa_external <- CESAnalysis(refset = "ces.refset.hg19")


# Create cesa object that will include only data in the paper
cesa_paper <- CESAnalysis(refset = "ces.refset.hg19")

# # Create cesa object that will include all data
cesa_all <- CESAnalysis(refset = "ces.refset.hg19")


# Check for duplicate data
maf_list <- list(dt1 = george_et_al_maf, dt2 = rudin_et_al_maf, dt3 = jiang_et_al_maf, dt4 = zhou_et_al_maf,
                 dt5 = song_et_al_maf, dt6 = chen_et_al_maf, dt7 = wang_et_al_maf)
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

cesa_external <- load_maf(cesa = cesa_external, maf = song_et_al_maf, maf_name = "song_et_al")
cesa_external <- load_sample_data(cesa_external, song_et_al_clinical)

cesa_all <- load_maf(cesa = cesa_all, maf = song_et_al_maf, maf_name = "song_et_al")
cesa_all <- load_sample_data(cesa_all, song_et_al_clinical)

cesa_external <- load_maf(cesa = cesa_external, maf = chen_et_al_maf, maf_name = "chen_et_al")
cesa_external <- load_sample_data(cesa_external, chen_et_al_clinical)

cesa_all <- load_maf(cesa = cesa_all, maf = chen_et_al_maf, maf_name = "chen_et_al")
cesa_all <- load_sample_data(cesa_all, chen_et_al_clinical)

cesa_external <- load_maf(cesa = cesa_external, maf = wang_et_al_maf, maf_name = "wang_et_al")
cesa_all <- load_maf(cesa = cesa_all, maf = wang_et_al_maf, maf_name = "wang_et_al")

# Summarize the amount of clinical and genomic data loaded in analysis

data_summary <- summarize_maf_sources(cesa_all)


save_cesa(cesa_external, "cesa_external.rds")
# cesa_external <- readRDS("cesa_external.rds")
save_cesa(cesa_paper, "cesa_paper_new.rds")
# cesa_paper <- readRDS("cesa_paper_new.rds")
save_cesa(cesa_all, "cesa_all.rds")
# cesa_all <- readRDS("cesa_all.rds")

# ---------------Assign COSMIC Mutational Signatures---------------------

# Define signature exclusions and use reference-based assignment (either "MutationalPatterns" or "deconstructSigs") to assign mutational signatures
assign_mutational_signatures("MutationalPatterns")

biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended==FALSE]
biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended==FALSE]

# Normalize TMB to mutations/Mb (all cohorts are WES; assume ~38 Mb exome)
exome_mb <- 38

biological_weights_paper[, TMB_per_Mb := total_snvs / exome_mb]
biological_weights_external[, TMB_per_Mb := total_snvs / exome_mb]
# 
# # non-bootstrapped optimal cutoff for paper data
# 
# # Load clinical data and prepare for analysis
# clinical_data_paper <- cesa_paper$samples %>%
#   filter(
#     !is.na(`Status.(at.time.of.last.follow-up)`),
#     !is.na(`overall_survival.(months)`)
#   ) %>%
#   mutate(
#     status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)
#   ) %>%
#   select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
# 
# # Merge mutation and clinical data
# paper_optimal_cutoff <- biological_weights_paper %>%
#   inner_join(clinical_data_paper, by = "Unique_Patient_Identifier")
# 
# # Perform optimal cutoff tests (non-bootstrap)
# SBS4_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
#   Surv(`overall_survival.(months)`, status_numeric) ~ SBS4,
#   data = paper_optimal_cutoff,
#   smethod = "LogRank", pmethod = "HL"
# )
# 
# SBS4_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
#   Surv(TMB_per_Mb) ~ SBS4,
#   data = paper_optimal_cutoff,
#   smethod = "LogRank", pmethod = "HL"
# )
# 
# SBS13_paper_optimal_survival_cutoff_for_paper_non_bs <- maxstat.test(
#   Surv(`overall_survival.(months)`, status_numeric) ~ SBS13,
#   data = paper_optimal_cutoff,
#   smethod = "LogRank", pmethod = "HL"
# )
# 
# SBS13_paper_optimal_TMB_cutoff_for_paper_non_bs <- maxstat.test(
#   Surv(TMB_per_Mb) ~ SBS13,
#   data = paper_optimal_cutoff,
#   smethod = "LogRank", pmethod = "HL"
# )

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
# Complete the following analyses: Kaplan-Meier analysis using 0.25 binary
SBS4_results <- analyze_signature_complete("SBS4")

# Chen-et-al did not have any samples with >0.25 SBS13 activity
# jiang_et_al did not have any samples with >0.25 SBS13 activity
SBS13_results <- analyze_signature_complete("SBS13")


# Hazard ratios w/ 95% CI

# SBS4
summary(SBS4_results$all$hazard_0.25) # 0.9976 [0.7487 - 1.329], n = 362
summary(SBS4_results$chen_et_al$hazard_0.25) # 0.606 [0.1443 - 2.544], n = 18
summary(SBS4_results$george_et_al$hazard_0.25)# 1.566 [0.9456 - 2.595], n = 101
summary(SBS4_results$jiang_et_al$hazard_0.25)# 1.143 [0.5519 - 2.367], n = 95
summary(SBS4_results$liu_et_al$hazard_0.25)# 0.5689 [0.3352 - 0.9655], n = 108
summary(SBS4_results$zhou_et_al$hazard_0.25) # 0.9736 [0.3458 - 2.741], n = 40

# SBS13
summary(SBS13_results$all$hazard_0.25) # 0.5143 [0.127 - 2.082], n = 362
# Chen-et-al did not have any samples with >0.25 SBS13 activity
summary(SBS13_results$george_et_al$hazard_0.25)# 1.066e-07 [0 - Inf], n = 101
# jiang_et_al did not have any samples with >0.25 SBS13 activity
summary(SBS13_results$liu_et_al$hazard_0.25)# 1.004 [0.1384 - 7.28], n = 108
summary(SBS13_results$zhou_et_al$hazard_0.25) # 1.202 [0.1525 - 9.47], n = 40

# Tried just using Li et al., 2025 signature data but could not replicate the survival findings
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

# Subset of paper_signatures in signature_analysis_paper

paper_signatures <- paper_signatures[paper_signatures$Unique_Patient_Identifier %in% survival_data_paper$Unique_Patient_Identifier,]
paper_sig_survival_data_paper <- survival_data_paper[survival_data_paper$Unique_Patient_Identifier %in% paper_signatures$Unique_Patient_Identifier, ]

# Merge clinical and signature data
# paper_signatures_optimal_cutoff <- paper_signatures %>%
#   inner_join(paper_sig_survival_data_paper, by = "Unique_Patient_Identifier") %>%
#   select(Unique_Patient_Identifier, SBS4, SBS5, SBS13, overall_survival_months, status_numeric)

# # Compute optimal cutoff using maxstat.test
# SBS4_paper_signatures_optimal_cutoff <- maxstat.test(
#   Surv(overall_survival_months, status_numeric) ~ SBS4, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
#   pmethod = "HL"
# )
# SBS13_paper_signatures_optimal_cutoff <- maxstat.test(
#   Surv(overall_survival_months, status_numeric) ~ SBS13, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
#   pmethod = "HL"
# )
# 
# SBS5_paper_signatures_optimal_cutoff <- maxstat.test(
#   Surv(overall_survival_months, status_numeric) ~ SBS5, data = paper_signatures_optimal_cutoff, smethod = "LogRank", 
#   pmethod = "HL"
# )

sig_result_SBS4 <- analyze_signature(signature_name = "SBS4", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)
summary(sig_result_SBS4$hazard_0.25)

sig_result_SBS13 <- analyze_signature(signature_name = "SBS13", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)
summary(sig_result_SBS13$hazard_0.25)

sig_result_SBS5 <- analyze_signature(signature_name = "SBS5", dataset = "paper signature assignment, paper", biological_weight_table = paper_signatures, survival_data = paper_sig_survival_data_paper)

# # ---------- Multivariable Cox regression per cohort --------------------
# 
# george_samples <- as.data.frame(cesa_paper$samples) %>%
#   
#   filter(maf_source == "george_et_al") %>%
#   
#   mutate(
#     age            = as.numeric(as.character(age)),
#     
#     OS_months      = as.numeric(
#       as.character(`overall_survival.(months)`)
#     ),
#     
#     status_numeric = ifelse(
#       `Status.(at.time.of.last.follow-up)` == "dead",
#       1,
#       0
#     ),
#     
#     prior_tx = as.integer(
#       previous.therapeutic.treatment.for.SCLC
#     ),
#     
#     sex = as.integer(sex == "male"),
#     
#     smoking = as.integer(
#       smoking_status == "smoker"
#     ),
#     
#     stage = case_when(
#       grepl("^I[Vv]", staging) ~ 4L,
#       grepl("^III", staging)  ~ 3L,
#       grepl("^II", staging)   ~ 2L,
#       grepl("^I", staging)    ~ 1L,
#       TRUE ~ NA_integer_
#     )
#   ) %>%
#   
#   filter(
#     !is.na(OS_months),
#     !is.na(status_numeric)
#   )
# 
# 
# george_cox_df <- biological_weights_paper %>%
#   
#   select(
#     Unique_Patient_Identifier,
#     SBS4,
#     SBS13
#   ) %>%
#   
#   inner_join(
#     george_samples,
#     by = "Unique_Patient_Identifier"
#   ) %>%
#   
#   mutate(
#     SBS4_high = as.integer(
#       SBS4 >= SBS4_survival_optimal_cutoff_paper_mean
#     ),
#     
#     SBS13_high = as.integer(
#       SBS13 >= SBS13_survival_optimal_cutoff_paper_mean
#     )
#   )
# 
# 
# # ---- George univariate models -----------------------------------------
# 
# cox_george_SBS4_uni <- run_cox(
#   george_cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS4_high",
#   "George_SBS4_uni"
# )
# 
# cox_george_SBS13_uni <- run_cox(
#   george_cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS13_high",
#   "George_SBS13_uni"
# )
# 
# 
# # ---- George multivariate models ---------------------------------------
# 
# cox_george_SBS4 <- run_cox(
#   george_cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS4_high + age + sex + stage + smoking + prior_tx",
#   "George_SBS4"
# )
# 
# cox_george_SBS13 <- run_cox(
#   george_cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS13_high + age + sex + stage + smoking + prior_tx",
#   "George_SBS13"
# )
# 
# 
# # ---------- External cohort --------------------------------------------
# 
# external_samples <- as.data.frame(cesa_external$samples) %>%
#   
#   mutate(
#     age = as.numeric(as.character(age)),
#     
#     OS_months = as.numeric(
#       as.character(`overall_survival.(months)`)
#     ),
#     
#     status_numeric = ifelse(
#       `Status.(at.time.of.last.follow-up)` == "dead",
#       1,
#       0
#     ),
#     
#     prior_tx = as.integer(
#       previous.therapeutic.treatment.for.SCLC
#     ),
#     
#     sex = as.integer(sex == "male"),
#     
#     smoking = as.integer(
#       smoking_status == "smoker"
#     ),
#     
#     stage = case_when(
#       grepl("^I[Vv]", staging) ~ 4L,
#       grepl("^III", staging)  ~ 3L,
#       grepl("^II", staging)   ~ 2L,
#       grepl("^I", staging)    ~ 1L,
#       TRUE ~ NA_integer_
#     )
#   ) %>%
#   
#   filter(
#     !is.na(OS_months),
#     !is.na(status_numeric)
#   )
# 
# 
# cox_df <- biological_weights_external %>%
#   
#   select(
#     Unique_Patient_Identifier,
#     SBS4,
#     SBS13
#   ) %>%
#   
#   inner_join(
#     all_samples,
#     by = "Unique_Patient_Identifier"
#   ) %>%
#   
#   mutate(
#     SBS4_high = as.integer(
#       SBS4 >= SBS4_survival_optimal_cutoff_paper_mean
#     ),
#     
#     SBS13_high = as.integer(
#       SBS13 >= SBS13_survival_optimal_cutoff_paper_mean
#     )
#   )
# 
# 
# # ---- Combined univariate models ---------------------------------------
# 
# cox_SBS4_uni <- run_cox(
#   cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS4_high",
#   "cox_SBS4_uni"
# )
# 
# cox_SBS13_uni <- run_cox(
#   cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS13_high",
#   "cox_SBS13_uni"
# )
# 
# 
# # ---- Combined multivariate models -------------------------------------
# 
# cox_SBS4 <- run_cox(
#   cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS4_high + age + stage + smoking",
#   "SBS4"
# )
# 
# cox_SBS13 <- run_cox(
#   cox_df,
#   "Surv(OS_months, status_numeric) ~ SBS13_high + age + stage + smoking",
#   "SBS13"
# )
# 
# # ---------- Multivariable Cox forest plots --------------------------------
# 
# make_forest_plot <- function(uni_cox, multi_cox, sig_name,
#                              x_limits = c(0.05, 20)) {
#   
#   font <- 11
#   
#   # automatic log-scale "nice" breaks:
#   # 1, 2, 5 × 10^k
#   all_breaks <- c(1, 2, 5) %>%
#     outer(10^(-10:10), `*`) %>%
#     as.vector() %>%
#     sort()
#   
#   hr_breaks <- all_breaks[
#     all_breaks >= x_limits[1] &
#       all_breaks <= x_limits[2]
#   ]
#   
#   # scientific notation labels parsed by ggplot
#   sci_labels <- sapply(hr_breaks, function(x) {
#     
#     exp <- floor(log10(x))
#     mant <- signif(x / (10^exp), 1)
#     
#     if (mant == 1) {
#       paste0("10^", exp)
#     } else {
#       paste0(mant, " %*% 10^", exp)
#     }
#   })
#   
#   var_labels <- c(
#     SBS4_high  = "High SBS4 activity",
#     SBS13_high = "High SBS13 activity",
#     age        = "Age",
#     sex        = "Sex (male)",
#     stage      = "Stage",
#     smoking    = "Smoking",
#     prior_tx   = "Prior treatment"
#   )
#   
#   extract_cox_table <- function(cox_result, model_label) {
#     
#     ci <- cox_result$fit_summary$conf.int
#     pv <- cox_result$fit_summary$coefficients[, "Pr(>|z|)"]
#     
#     data.frame(
#       Variable     = rownames(ci),
#       Hazard_Ratio = ci[, "exp(coef)"],
#       CI_lower     = ci[, "lower .95"],
#       CI_upper     = ci[, "upper .95"],
#       p_value      = pv,
#       Model        = model_label,
#       row.names    = NULL,
#       stringsAsFactors = FALSE
#     )
#   }
#   
#   combined_table <- bind_rows(
#     extract_cox_table(uni_cox,   "Univariate"),
#     extract_cox_table(multi_cox, "Multivariate")
#   ) %>%
#     mutate(
#       Variable = dplyr::recode(Variable, !!!var_labels),
#       Model = factor(
#         Model,
#         levels = c("Univariate", "Multivariate")
#       )
#     )
#   
#   var_order <- combined_table %>%
#     mutate(is_sig = grepl(sig_name, Variable)) %>%
#     arrange(desc(is_sig), Variable) %>%
#     pull(Variable) %>%
#     unique()
#   
#   combined_table$Variable <- factor(
#     combined_table$Variable,
#     levels = rev(var_order)
#   )
#   
#   combined_table <- combined_table %>%
#     mutate(
#       CI_lower_plot = pmax(CI_lower, x_limits[1]),
#       CI_upper_plot = pmin(CI_upper, x_limits[2]),
#       
#       lower_truncated = CI_lower < x_limits[1],
#       upper_truncated = CI_upper > x_limits[2]
#     )
#   
#   ggplot(
#     combined_table,
#     aes(x = Hazard_Ratio, y = Variable, color = Model)
#   ) +
#     
#     geom_point(
#       position = position_dodge(width = 0.5),
#       size = font * 0.18
#     ) +
#     
#     geom_errorbarh(
#       aes(
#         xmin = CI_lower_plot,
#         xmax = CI_upper_plot
#       ),
#       position = position_dodge(width = 0.5),
#       height = 0.15,
#       linewidth = 0.5
#     ) +
#     
#     geom_vline(
#       xintercept = 1,
#       linetype = "dashed",
#       linewidth = 0.5
#     ) +
#     
#     scale_x_log10(
#       limits = x_limits,
#       breaks = hr_breaks,
#       labels = parse(text = sci_labels),
#       expand = expansion(mult = c(0.02, 0.05))
#     ) +
#     
#     labs(
#       x = "Hazard Ratio (log scale)",
#       y = NULL,
#       title = paste0(
#         sig_name,
#         " — univariate vs. multivariate Cox (n=",
#         multi_cox$n,
#         ")"
#       )
#     ) +
#     
#     theme_minimal(base_size = font) +
#     
#     theme(
#       text = element_text(size = font),
#       
#       plot.title = element_text(
#         hjust = 0.5,
#         size = font * 1.1,
#         margin = margin(b = 4)
#       ),
#       
#       axis.title.x = element_text(
#         margin = margin(t = 8)
#       ),
#       
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = font * 0.85
#       ),
#       
#       axis.text.y = element_text(
#         size = font * 0.9,
#         margin = margin(r = 4)
#       ),
#       
#       legend.position = "top",
#       legend.title = element_blank(),
#       
#       legend.text = element_text(
#         size = font * 0.9
#       ),
#       
#       legend.margin = margin(t = -6, b = -8),
#       legend.box.margin = margin(0, 0, 0, 0),
#       legend.spacing.x = unit(4, "pt"),
#       
#       panel.grid.minor = element_blank(),
#       panel.grid.major.y = element_blank(),
#       
#       plot.margin = margin(
#         t = 4,
#         r = 4,
#         b = 20,
#         l = 4
#       )
#     ) +
#     
#     guides(
#       color = guide_legend(nrow = 1)
#     )
# }
# 
# forest_SBS4_george <- make_forest_plot(
#   cox_george_SBS4_uni,
#   cox_george_SBS4,
#   "SBS4"
# )
# 
# forest_SBS4 <- make_forest_plot(
#   cox_SBS4_uni,
#   cox_SBS4,
#   "SBS4"
# )
# 
# forest_SBS13_george <- make_forest_plot(
#   cox_george_SBS13_uni,
#   cox_george_SBS13,
#   "SBS13",
#   x_limits = c(1e-06, 10)
# )
# 
# forest_SBS13 <- make_forest_plot(
#   cox_SBS13_uni,
#   cox_SBS13,
#   "SBS13"
# )
# 
# 
# # george_ethnicity_df <- george_cox_df %>%
# #   filter(ethnicity %in% c("asian", "caucasian"), !is.na(OS_months), !is.na(status_numeric))
# # 
# # # Interaction model: tests whether SBS4-survival slope differs by ethnicity
# # sbs4_ethnicity_interaction <- coxph(
# #   Surv(OS_months, status_numeric) ~ SBS4 * ethnicity ,
# #   data = george_ethnicity_df
# # )
# # sbs4_ethnicity_interaction_p <- summary(sbs4_ethnicity_interaction)$coefficients["SBS4:ethnicitycaucasian", "Pr(>|z|)"]
# # 
# # # Stratified univariate Cox: SBS4 HR within each ethnicity group
# # george_asian      <- george_ethnicity_df %>% filter(ethnicity == "asian")
# # george_caucasian  <- george_ethnicity_df %>% filter(ethnicity == "caucasian")
# # 
# # cox_sbs4_asian     <- run_cox(george_asian,     "Surv(OS_months, status_numeric) ~ SBS4", "Asian")
# # cox_sbs4_caucasian <- run_cox(george_caucasian, "Surv(OS_months, status_numeric) ~ SBS4", "Caucasian")
# # 
# # # Forest plot: SBS4 HR by ethnicity
# # sbs4_ethnicity_table <- bind_rows(
# #   extract_cox_table(cox_sbs4_asian,     "Asian"),
# #   extract_cox_table(cox_sbs4_caucasian, "Caucasian")
# # ) %>%
# #   mutate(Model = factor(Model, levels = c("Asian", "Caucasian")))
# # 
# # forest_SBS4_ethnicity <- ggplot(sbs4_ethnicity_table,
# #                                 aes(x = Hazard_Ratio, y = Model, color = Model)) +
# #   geom_point(size = 11 * 0.25) +
# #   geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15, linewidth = 0.5) +
# #   geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
# #   scale_x_log10(expand = expansion(mult = c(0.05, 0.05))) +
# #   annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
# #            label = paste0("Interaction P = ", signif(sbs4_ethnicity_interaction_p, 2)),
# #            size = 11 / 3, color = "grey30") +
# #   labs(x = "HR per unit SBS4 activity (log scale)", y = NULL,
# #        title = "SBS4 vs. OS by ethnicity (George et al.)") +
# #   theme_minimal(base_size = 11) +
# #   theme(legend.position = "none",
# #         panel.grid.minor   = element_blank(),
# #         panel.grid.major.y = element_blank(),
# #         plot.margin        = margin(t = 4, r = 8, b = 8, l = 4))

# ---------- Annotation discrepancy ------
# There is a difference between our attempted replication (sig_result_SBS4$plot_0.25 and sig_result_SBS13$plot_0.25) with the published Fig 2D and 2F (Li et al., 2025). Goal is to identify which annotations were modified in the published results.

# George et al. vital status
george_ground_truth <- george_et_al_clinical %>%
  select(Unique_Patient_Identifier, `overall_survival.(months)`, `Status.(at.time.of.last.follow-up)`) %>% filter(!is.na(`Status.(at.time.of.last.follow-up)`))

# Li et al. signature activities and group assignments (Table S3)
li_groups <- paper_signatures %>%
  select(Unique_Patient_Identifier, li_SBS4 = SBS4, li_SBS13 = SBS13) %>%
  mutate(
    li_SBS4_group  = ifelse(li_SBS4  >= 0.25, "high", "low"),
    li_SBS13_group = ifelse(li_SBS13 >= 0.25, "high", "low")
  )

discrepancy_tbl <- li_groups %>%
  inner_join(george_ground_truth, by = "Unique_Patient_Identifier") %>%
  arrange(`overall_survival.(months)`)

# Flag patients at suspicious OS time points determined by visually comparing the attempted replicated figures with Li et al., 2025 figures
suspicious_SBS13_OS <- c(23, 30, 46)
suspicious_SBS4_OS  <- c(0, 20, 23, 72)

discrepancy_tbl <- discrepancy_tbl %>%
  mutate(
    flag_SBS13 = li_SBS13_group == "high" &
      sapply(`overall_survival.(months)`, function(t) any(abs(t - suspicious_SBS13_OS) <= 0)),
    flag_SBS4  = li_SBS4_group  == "low"  &
      sapply(`overall_survival.(months)`, function(t) any(abs(t - suspicious_SBS4_OS)  <= 0))
  )

cat("\n=== SBS13-high patients at OS ≈ 23/30/46 months (George et al. status) ===\n")
SBS13_discrepancies <- discrepancy_tbl %>%
    filter(flag_SBS13) %>%
    select(
      Unique_Patient_Identifier,
      li_SBS13_group,
      `overall_survival.(months)`,
      `Status.(at.time.of.last.follow-up)`
    ) %>%
    mutate(
      li_status = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", "alive", "dead")
    ) %>%
    group_by(`overall_survival.(months)`) %>%
    summarise(
      potential_pt_ids = paste(
        Unique_Patient_Identifier,
        collapse = " or "
      ),
      li_SBS13_group = dplyr::first(li_SBS13_group),
      li_SBS4_group = NA,
      george_status = dplyr::first(`Status.(at.time.of.last.follow-up)`),
      li_status = dplyr::first(li_status),
      .groups = "drop"
    )


cat("\n=== SBS4-low patients at OS ≈ 0/20/23/72 months (George et al. status) ===\n")
SBS4_discrepancies <- discrepancy_tbl %>%
        filter(flag_SBS4) %>%
        select(
          Unique_Patient_Identifier,
          li_SBS4_group,
          `overall_survival.(months)`,
          `Status.(at.time.of.last.follow-up)`
        ) %>%
        mutate(
          li_status = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", "alive", "dead")
        ) %>% filter(`Status.(at.time.of.last.follow-up)` == "dead") %>%
        group_by(`overall_survival.(months)`) %>%
        summarise(
          potential_pt_ids = paste(
            Unique_Patient_Identifier,
            collapse = " or "
          ),
          li_SBS4_group = dplyr::first(li_SBS4_group),
          li_SBS13_group = NA,
          george_status = dplyr::first(`Status.(at.time.of.last.follow-up)`),
          li_status = dplyr::first(li_status),
          .groups = "drop"
        )

write.csv(rbind(SBS4_discrepancies, SBS13_discrepancies), "annotation_discrepancy_table.csv", row.names = FALSE)
